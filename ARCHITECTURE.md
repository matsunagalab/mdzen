# MCP-MD アーキテクチャ・実装プラン

## 1. プロジェクト概要

### 目的とポジショニング

**Amber系に最適化したAIエージェント＋MCPツール群**

- **主軸**: Amber/GAFF/OpenFF/ParmEd/OpenMM エコシステムに特化
- **非競合**: CHARMM-GUIとは棲み分け（CHARMM系は変換経由で二次対応、将来拡張）
- **永続化**: MCP標準でツール接続を維持可能（将来のLLM/実行基盤の更新に強い）
- **ホスト/クライアント**: [Strands Agents](https://github.com/Strands-AI/strands)に統一（MCPクライアント統合）

### 主要技術スタック

- **Strands Agents**: 永続AIエージェントフレームワーク（MCPクライアント内蔵）
- **Boltz-2**: 構造予測・複合体生成ツール
- **AmberTools**: 完全OSS、配位子パラメータ化（GAFF2 + AM1-BCC）
- **OpenMM**: Pythonプログラマブル、GPU最適化、プロダクション対応MD
- **MCP (Model Context Protocol)**: 標準化されたツール統合（ツールの永続性・相互運用性）

### 主要機能

1. **Hybrid Planning**: 固定スケルトン + 自律サブルーチン（意思決定ログ記録）
2. **品質保証**: 自作MolProbity等による物理化学的一貫性チェック
3. **再現性**: Plan/決定/生成物をJSON保存

---

## 2. 全体アーキテクチャ（ハイブリッド設計）

### システム構成

```
[Chat UI / Jupyter / CLI]
    ↓
[Strands Agent] ──────────┐
  ├─ Planner              │ (永続メモリ)
  │   └─ 固定スケルトン    │  - ユーザ既定（pH, 塩, box）
  ├─ Memory               │  - 過去の実行履歴
  │   └─ User Preferences │  - 決定根拠ログ
  ├─ Policy               │
  │   └─ 自律サブルーチン   │
  └─ FastMCP Client ──────├─── [FastMCP Servers] 🆕
                          │
                          ├─ Structure Server
                          │   ├─ fetch_pdb
                          │   ├─ clean_structure
                          │   ├─ add_hydrogens
                          │   ├─ protonate_structure
                          │   ├─ detect_modifications
                          │   └─ validate_structure
                          │
                          ├─ Genesis Server 🆕
                          │   ├─ boltz2_protein_from_seq
                          │   ├─ boltz2_protein_from_fasta
                          │   └─ boltz2_multimer
                          │
                          ├─ Complex Server 🆕
                          │   ├─ boltz2_complex
                          │   ├─ boltz2_screen_ligands
                          │   ├─ smina_dock
                          │   └─ refine_poses
                          │
                          ├─ Ligand Server
                          │   ├─ smiles_to_3d
                          │   ├─ generate_gaff_params
                          │   ├─ create_ligand_lib
                          │   └─ parameterize_ligand_complete
                          │
                          ├─ Assembly Server
                          │   ├─ build_system_tleap
                          │   ├─ build_membrane_system
                          │   └─ build_mixed_solvent
                          │
                          ├─ Export Server
                          │   ├─ export_amber
                          │   ├─ export_gromacs
                          │   ├─ export_openmm
                          │   ├─ package_system
                          │   └─ convert_format
                          │
                          └─ QC/Min Server 🆕
                              ├─ openmm_minimize
                              ├─ clash_check
                              ├─ bond_check
                              ├─ chirality_check
                              ├─ run_full_qc
                              └─ posebusters_check

    ↓
[Persistent Storage]
  └─ runs/<timestamp>/
      ├─ plan.json        (固定スケルトン + 決定ログ)
      ├─ outputs/         (PDB, prmtop, inpcrd, etc.)
      ├─ qc_report.json   (PoseBusters, 最小化指標)
      └─ metadata.json    (seed, hash, 再現用)
```

### FastMCP統合の特徴

- **モジュラー設計**: 各サーバーファイルが完全に独立して動作可能
- **自動スキーマ生成**: 型ヒントとdocstringから自動的にMCPツールスキーマを生成
- **標準準拠**: MCP標準プロトコルに完全準拠
- **開発効率**: デコレータベースのシンプルなAPI（`@mcp.tool`）
- **独立実行**: 各サーバーが `python -m servers.{server_name}` で単独起動可能
- **共通ライブラリ**: `common/` モジュールで外部ツール実行とユーティリティを共有

### ハイブリッド設計の核心

#### 1. 固定スケルトン（Plan）
```
fetch/generate → repair/protonate → ligand_param → 
complex → assemble → solvate/ions → export → 
minimize → package
```

**特徴**:
- 工程は固定（Amber特化パイプライン）
- 順序は決定論的
- 学術的再現性を担保

#### 2. 自律サブルーチン（Policy）
各工程内でツール・パラメータを動的選択：

```python
# 例: 複合体生成の意思決定
if pdb_exists:
    tool = "fetch_pdb"
elif fasta_provided:
    tool = "boltz2_protein_from_seq"
    log_decision("Genesis from sequence", reason="No PDB available")

# 複合体ポーズ生成
if use_ai_model:
    poses = boltz2_complex(protein, ligand, top_k=5)
    log_decision("Boltz-2 complex", affinity=poses[0].affinity)
    
    if refine_poses:
        poses = smina_dock(poses, local_search=True)
        log_decision("Smina refinement", reason="Improve local geometry")
```

**特徴**:
- 失敗時は自動で1回再試行
- それでもNGなら3択程度の簡潔な質問
- すべての決定をJSON記録

---

## 3. MCPツール統合

### Genesis MCP: 構造生成

FASTA配列からPDB構造を生成：

```python
# FASTA → PDB
protein_pdb = boltz2_protein_from_seq(
    sequence="MKTAYIAKQRQISFVKSHFSRQ...",
    num_models=5
)
```

### Complex MCP: 複合体生成

タンパク質-配位子複合体の姿勢予測：

```python
# 受容体 + SMILES → 複合体候補
complexes = boltz2_complex(
    protein_pdb="receptor.pdb",
    ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
    top_k=10
)

# Sminaで局所精密化（オプション）
refined_poses = smina_dock(
    receptor="receptor.pdb",
    ligands=complexes[:5],
    local_search=True
)
```

### QC/Min MCP: 品質保証

物理化学的一貫性チェック：

```python
# PoseBustersチェック
qc_report = posebusters_check(pdb_file="complex.pdb")

# OpenMM最小化
minimized = openmm_minimize(
    prmtop="system.prmtop",
    inpcrd="system.inpcrd",
    max_iterations=5000
)
```

---

## 4. 永続エージェント（Strands × MCP）

### Strands Agentsとは

- **公式**: https://github.com/Strands-AI/strands
- **特徴**: 永続メモリ、マルチターン対話、MCP統合、セッション管理
- **MCP統合**: `with mcp_client:` コンテキストで安全にツール呼び出し

### 運用のキーポイント

#### 1. MCPセッション管理
```python
from strands import Agent
from strands.mcp import MCPClient

agent = Agent(
    name="md-assistant",
    model="gpt-4o",  # またはローカルLLM
    mcp_clients=[
        MCPClient("structure_server"),
        MCPClient("complex_server"),
        MCPClient("assembly_server"),
        # ...
    ]
)

async with agent:
    result = await agent.execute(
        "Fetch PDB 1ABC and dock Aspirin"
    )
```

#### 2. 永続メモリ
```python
# ユーザ既定の保存
agent.memory.set("user_preferences", {
    "ph": 7.4,
    "salt_concentration": 0.15,  # M
    "water_model": "TIP3P",
    "force_field": "ff19SB",
    "known_binding_sites": ["SER195", "HIS57", "ASP102"]
})

# 過去の実行履歴
agent.memory.add("execution_history", {
    "timestamp": "2025-01-20T10:30:00Z",
    "query": "PDB 1ABC + Aspirin",
    "success": True,
    "output": "runs/20250120_103000/"
})
```

#### 3. 意思決定の記録
```python
# すべての決定をログ
def log_decision(step: str, tool: str, params: dict, reason: str):
    agent.memory.append("decisions", {
        "timestamp": datetime.utcnow().isoformat(),
        "step": step,
        "tool": tool,
        "params": params,
        "reason": reason
    })

# 例
log_decision(
    step="complex_generation",
    tool="boltz2_complex",
    params={"top_k": 5, "use_msa": True},
    reason="High confidence structure needed, MSA available"
)
```

#### 4. 安全性・認証
- **ID分散の課題**: MCP標準だが認証・権限は別管理が推奨
- **実装方針**:
  - Strands側で機密値（API key等）を秘匿
  - 一時クレデンシャル運用
  - ローカル実行前提（外部API最小化）

---

## 5. ユーザ体験（UX）設計

### 典型的な対話フロー

#### 1. 最初の返答
```
User: "PDB 1ABCにAspirinをドッキングしてMD系を作って"

Agent:
了解しました。以下のプランで進めます：

【固定スケルトン】
1. PDB取得（1ABC）
2. 構造修復・プロトネーション
3. リガンド処理（Aspirin, GAFF2/AM1-BCC）
4. 複合体生成（Boltz-2 + Smina）
5. 系構築（tleap）
6. 溶媒化・イオン付与
7. Amber出力（prmtop/inpcrd）
8. 最小化・QC
9. パッケージング

【質問: 以下を教えてください】
- pH: [default: 7.4]
- 塩濃度: [default: 0.15 M]
- Box padding: [default: 12.0 Å]
- 出力形式: [Amber / GROMACS / OpenMM]
- 既知結合部位: [あれば指定、なければBoltz-2で推定]
```

#### 2. 実行中の可視化
```
[Step 1/9] PDB取得 ✅ (1ABC.pdb, 1234 atoms)
[Step 2/9] 構造修復 ✅ (欠損残基 3箇所補完)
[Step 3/9] リガンド処理 ⏳
  - SMILES → 3D (RDKit) ✅
  - GAFF2パラメータ化 (AM1-BCC) ⏳
    Decision: AM1-BCC選択（バランス重視、計算時間 < 1min）

[中間プレビュー]
🔗 http://localhost:8080/view/intermediate.pdb (NGLビューワ)
```

#### 3. 失敗時の誘導
```
[Step 4/9] 複合体生成 ⚠️ エラー
  - Boltz-2予測: 親和性が極めて低い（binder_prob = 0.12）
  
【自動リトライ】
  - Smina局所サーチで代替候補を探索中...
  - 結果: 候補なし

【提案: 以下から選択してください】
a) 結合サイトを手動指定（推奨残基: SER195, HIS57, ASP102）
b) リガンドの3Dコンフォーマを変更（ETKDG → UFF）
c) Boltz-2の設定変更（MSA使用、top_k=10）
```

---

## 6. ロードマップ（Amber特化 → 拡張）

### Phase 1: MVP（最小実装）🎯

**目標**: Amber特化の基本パイプライン

**機能**:
- PDB/FASTA → Boltz-2 or fetch → 修復 → GAFF2/AM1-BCC
- Boltz-2複合体 ± Smina → tleap → Amber出力（prmtop/inpcrd）
- 基本的なQC（最小化、電荷整合）

**MCPサーバー**:
- Structure MCP
- Genesis MCP
- Complex MCP
- Ligand MCP
- Assembly MCP
- Export MCP（Amberのみ）

**期間**: 4-6週間

### Phase 2: QC強化

**目標**: 学術発表レベルの品質保証

**追加機能**:
- 自作MolProbityによる物理化学的一貫性チェック
- QCレポート生成（JSON + Markdown）

**期間**: 2-3週間

### Phase 3: 出力拡張

**目標**: GROMACS/OpenMM対応、膜系

**追加機能**:
- GROMACS出力（ParmEd）
- OpenMM XML出力
- 膜系構築（Packmol-Memgen）
- 混合溶媒対応（Packmol）

**期間**: 3-4週間

### Phase 4: HPC/永続運用

**目標**: 大規模スクリーニング、長期運用

**追加機能**:
- Strandsエージェントのキュー実行
- 結果キャッシュ（同一入力の再利用）
- 結果索引（検索可能なDB）
- HPC連携（Slurmジョブ投入）

**期間**: 4-6週間

### Phase 5: 将来拡張（低優先度）

**CHARMM系対応**:
- CHARMMパラメータ化（CGenFF）
- CHARMM-GUI出力への変換

**特殊系対応**:
- 糖鎖（GLYCAM）
- 金属中心（MCPB.py）
- RNA特化（OL3）

**プラグインアーキテクチャ**:
- 外部開発者がMCPサーバーを追加可能
- コミュニティ貢献の促進

---

## 7. 現在の実装状況（FastMCP統合完了）

### 実装済み（7 FastMCP Servers）✅

| Component | Status | 主要機能 |
|-----------|--------|---------|
| Structure Server | ✅ | PDB取得、PDBFixer、PDB2PQR、構造検証 |
| Genesis Server | ✅ 🆕 | Boltz-2タンパク質生成（FASTA→PDB、マルチマー） |
| Complex Server | ✅ 🆕 | Boltz-2複合体予測、Sminaドッキング、ポーズ精密化 |
| Ligand Server | ✅ | RDKit 3D生成、AmberTools GAFF2/AM1-BCC |
| Assembly Server | ✅ | tleap系構築、Packmol-Memgen膜系 |
| Export Server | ✅ | Amber/GROMACS/OpenMM形式変換、パッケージング |
| QC/Min Server | ✅ 🆕 | OpenMM最小化、衝突検出、結合長・キラリティチェック |

### FastMCP統合アーキテクチャ

#### 1. ディレクトリ構造
```
mcp-md/
├── common/              # 共通ライブラリ（新規）
│   ├── __init__.py
│   ├── base.py         # BaseToolWrapper
│   └── utils.py        # 共通ユーティリティ
├── servers/            # FastMCPサーバー（7ファイル）
│   ├── __init__.py
│   ├── structure_server.py
│   ├── genesis_server.py
│   ├── complex_server.py
│   ├── ligand_server.py
│   ├── assembly_server.py
│   ├── export_server.py
│   └── qc_min_server.py
├── core/               # エージェント実装
│   ├── strands_agent.py    # FastMCP Client統合
│   ├── workflow_skeleton.py
│   ├── decision_logger.py
│   └── models.py
└── pyproject.toml      # fastmcp>=0.1.0追加
```

#### 2. FastMCP統合の実装パターン

各サーバーは以下の標準パターンで実装：

```python
from fastmcp import FastMCP
from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)
mcp = FastMCP("Server Name")

# 外部ツールラッパー初期化
tool_wrapper = BaseToolWrapper("tool_name", conda_env="mcp-md")

@mcp.tool
def tool_name(param1: str, param2: int = 0) -> dict:
    """Tool description
    
    Args:
        param1: Parameter description
        param2: Optional parameter
    
    Returns:
        Result dictionary
    """
    # 実装コード
    return result

if __name__ == "__main__":
    mcp.run()  # STDIO transport (default)
```

#### 3. Strands Agent → FastMCP Client統合

```python
from fastmcp import Client

# MCP設定（標準フォーマット）
mcp_config = {
    "mcpServers": {
        "structure": {
            "command": "python",
            "args": ["-m", "servers.structure_server"]
        },
        # ... 他のサーバー
    }
}

# FastMCP Clientで全サーバーに接続
async with Client(mcp_config) as client:
    tools = await client.list_tools()
    agent = Agent(model=model, tools=tools)
    response = await agent(user_query)
```

### 削除されたファイル（FastMCPに置き換え）

- ~~`tools/`~~ ディレクトリ全体（10ファイル） → `common/`に統合
- ~~`servers/base_server.py`~~ → FastMCP標準機能で代替
- ~~`servers/archive/`~~ → 旧実装削除

---

## 8. 技術詳細（Phase別）

以下は既存実装の技術詳細です。新設計への移行時に参照してください。

### Phase 1: Structure Server

**実装ファイル**:
- `servers/structure_server.py` (494行)
- `tools/boltz2_wrapper.py` (210行)
- `tools/pdbfixer_wrapper.py` (127行)
- `tools/pdb2pqr_wrapper.py` (129行)

**主要ツール**:
1. `fetch_pdb`: PDB/AlphaFold/PDB-REDO取得
2. `predict_structure_boltz2`: FASTA→構造（→ Genesis MCPに移行）
3. `predict_complex_with_affinity`: FASTA+SMILES→複合体（→ Complex MCPに移行）
4. `clean_structure`: PDBFixerクリーニング
5. `protonate_structure`: PDB2PQR+PROPKA
6. `detect_modifications`: ジスルフィド・金属検出

**リファクタリング**:
- Boltz-2関連ツールを Genesis/Complex MCP に分離
- Structure MCP は構造取得・修復・プロトネーションに専念

### Phase 2: Ligand Server

**実装ファイル**:
- `servers/ligand_server.py` (187行)
- `tools/rdkit_wrapper.py` (88行)
- `tools/ambertools_wrapper.py` (223行)

**主要ツール**:
1. `smiles_to_3d`: SMILES → 3D（RDKit ETKDG）
2. `generate_gaff_params`: antechamber + parmchk2（GAFF2/AM1-BCC）
3. `create_ligand_lib`: tleap用ライブラリ

**新設計での位置づけ**:
- そのまま維持（Amber特化の核心部分）
- OpenFF統合は Phase 3以降

### Phase 3: Docking Server

**実装ファイル**:
- `servers/docking_server.py` (135行)
- `tools/smina_wrapper.py` (162行)

**主要ツール**:
1. `prepare_receptor/ligand`: PDBQT変換
2. `dock_ligand_smina`: Sminaドッキング
3. `align_to_reference`: 既知リガンド整列

**新設計での位置づけ**:
- Complex MCP に統合
- Boltz-2複合体予測の補助ツールとして位置づけ

### Phase 4: Assembly Server

**実装ファイル**:
- `servers/assembly_server.py` (156行)
- `tools/ambertools_wrapper.py` - tleap統合
- `tools/packmol_wrapper.py` (144行)

**主要ツール**:
1. `build_system_tleap`: 完全MD系構築
2. `build_membrane_system`: Packmol-Memgen膜系

**新設計での位置づけ**:
- そのまま維持（Amber特化の核心）

### Phase 5: Protocol Server

**実装ファイル**:
- `servers/protocol_server.py` (220行)
- `tools/openmm_wrapper.py` (125行)

**主要ツール**:
1. `generate_openmm_minimization`: 最小化
2. `generate_openmm_equilibration`: 平衡化
3. `generate_openmm_production`: プロダクションMD

**新設計での位置づけ**:
- そのまま維持
- 最小化機能は QC/Min MCP にも複製

### Phase 6: Export Server

**実装ファイル**:
- `servers/export_server.py` (178行)
- ParmEd統合

**主要ツール**:
1. `export_amber`: prmtop/inpcrd
2. `export_gromacs`: ParmEd変換
3. `export_openmm`: XML
4. `package_system`: ZIP化

**新設計での位置づけ**:
- そのまま維持
- Phase 1はAmberのみ、Phase 3でGROMACS/OpenMM追加

---

## 9. 外部ツール依存関係

### conda経由（推奨）

```bash
conda create -n mcp-md python=3.11
conda activate mcp-md
conda install -c conda-forge ambertools packmol smina pdbfixer
```

### pip経由

```bash
pip install -e .  # pyproject.toml参照

# 主要パッケージ:
# - boltz>=2.0.0
# - pdb2pqr>=3.1.0, propka>=3.5.0
# - rdkit>=2023.9.1
# - openmm>=8.3.1, parmed>=4.3.0
# - openai>=1.0.0 (LM Studio用)
# - strands>=0.1.0 (Strands Agents)
# - mcp>=1.18.0
```

### Strands Agents

```bash
pip install strands-ai
```

---

## 10. 参考資料

### 主要論文

#### Boltz-2
```bibtex
@article{passaro2025boltz2,
  author = {Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and others},
  title = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
  journal = {bioRxiv},
  year = {2025}
}
```

### 外部リンク

- **Boltz-2**: https://github.com/jwohlwend/boltz
- **Strands Agents**: https://github.com/Strands-AI/strands
- **AmberTools**: https://ambermd.org/AmberTools.php
- **OpenMM**: https://openmm.org/
- **PoseBusters**: https://github.com/maabuu/posebusters
- **MCP Protocol**: https://modelcontextprotocol.io/

---

## 11. まとめ

### 現在地

- ✅ 6つのMCPサーバー実装済み（基本機能）
- ✅ 9つのツールラッパー完成
- ❌ Strands統合未実装
- ❌ Genesis/Complex MCP未実装
- ❌ QC/Min MCP未実装

### 次のステップ

1. **Strands Agent統合**（最優先、2週間）
2. **Genesis/Complex MCP実装**（2週間）
3. **QC/Min MCP実装**（PoseBusters統合、1週間）
4. **MVP完成**（Phase 1完了、4週間）

### プロジェクト特性

- **非競合**: CHARMM-GUIと棲み分け（Amber特化）
- **将来性**: MCP標準でツール永続化、LLM/実行基盤の更新に強い
