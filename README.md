# MCP-MD: Molecular Dynamics Input File Generation Agent

CHARMM-GUIに代わる、お手軽でフレクシブルなMD入力ファイル生成システム。Boltz-2による構造・親和性予測、AmberToolsによる配位子パラメータ化、OpenMMによるMD実行を統合。

## 特徴

- **Boltz-2統合**: FASTAやSMILESから高精度な構造予測と結合親和性予測
  - 構造予測、複合体+親和性予測、バーチャルスクリーニング
- **AmberTools完結**: 配位子パラメータ化に外部QMソフト不要（AM1-BCC電荷計算）
  - SMILES → 3D → GAFF2パラメータ → tleapライブラリの完全自動化
- **高度な構造処理**:
  - PDB2PQR+PROPKAによるpH指定プロトネーション
  - ジスルフィド結合・金属サイト自動検出
- **膜タンパク質系**: Packmol-Memgen統合で脂質二重層自動構築
- **OpenMM専用**: Pythonプログラマブルなプロダクションレディなスクリプト生成
- **LangGraph統合** 🆕: ステートフルなワークフロー、永続化、人間フィードバック
  - LangChain 1.0準拠のStateGraphベースの実装
  - langchain-mcp-adaptersで公式MCP統合
  - チェックポイント機能で中断・再開可能
- **FastMCP統合**: モジュラーな独立サーバー、型安全な自動スキーマ生成
  - 7つの独立したFastMCPサーバー(各サーバーが単独で動作可能)
  - デコレータベースのシンプルなAPI（`@mcp.tool`）
  - 標準MCP準拠で将来のLLM/実行基盤更新に強い

## 📚 ドキュメント

- **[ARCHITECTURE.md](ARCHITECTURE.md)** - プロジェクト全体のアーキテクチャ・実装プラン・技術仕様
- **[AGENTS.md](AGENTS.md)** - Cursor AI Agent設定とガイドライン
- **[.cursor/rules/](.cursor/rules/)** - プロジェクトルールと開発ワークフロー
  - **[project-rules.md](.cursor/rules/project-rules.md)** - プロジェクト全体のルール
  - **[notebook-development.md](.cursor/rules/notebook-development.md)** - Notebook開発ルール

## セットアップ

### 前提条件

- Python 3.11以上
- [conda](https://docs.conda.io/en/latest/) または [mamba](https://mamba.readthedocs.io/) (推奨)
- [LM Studio](https://lmstudio.ai/) (ローカルLLM実行)
- GPU推奨（Boltz-2、OpenMM高速化）

### インストール手順

#### 1. リポジトリのクローン

```bash
git clone https://github.com/matsunagalab/mcp-md.git
cd mcp-md
```

#### 2. conda環境のセットアップ（推奨）

```bash
# 1. conda環境作成
conda create -n mcp-md python=3.11
conda activate mcp-md

# 2. 外部ツールをインストール（conda-forge）
conda install -c conda-forge ambertools packmol smina pdbfixer

# 3. Python依存関係をインストール
pip install -e .

# 4. Boltz-2インストール（GPU版）
pip install "boltz[cuda]" -U
```

> **注意**: MCPサーバーやスクリプトを実行する際は、必ず`conda activate mcp-md`で環境を有効化してください。

#### 3. LM Studioのセットアップ

1. [LM Studio](https://lmstudio.ai/)をダウンロード・インストール
2. LM Studio GUIでモデルをダウンロード（推奨: `gemma-3-12b`）
3. `Local Server`タブで`Start Server`をクリック（デフォルト: `http://localhost:1234`）
4. 環境変数を設定（オプション）:

```bash
# ~/.bashrc または ~/.zshrc に追加
export LM_STUDIO_BASE_URL="http://localhost:1234/v1"
export LM_STUDIO_MODEL="gemma-3-12b"
```

> **ヒント**: 環境変数を設定しない場合、デフォルト値が使用されます。

## 使用方法

### 🚀 対話型チャット（推奨）

LangGraphの対話型ワークフローを使用：

```bash
# conda環境をアクティベート
conda activate mcp-md

# LM Studioを起動（別ターミナル）
# http://localhost:1234 でサーバーが起動していることを確認

# 対話型チャットを開始
mcp-md chat

# または、モデルやURLを指定
mcp-md chat --model gemma-3-12b --lm-studio-url http://192.168.1.100:1234/v1
```

チャット内で自然言語でリクエストを送信：

```
> Generate a protein structure from this FASTA sequence: MKFLKFSLLTAVLLSVVFAFSSCGDDDDTYPYDVPDYAG

> Create an MD system for protein 1ABC with ligand CCO (ethanol)

> Quality check my PDB file: structure.pdb
```

すべての実行状態は `checkpoints/workflow.db` に永続化され、中断・再開が可能：

```
> resume <thread_id>     # 中断したワークフローを再開
> history <thread_id>    # ワークフローの実行履歴を表示
```

### MCP Inspectorでデバッグ（開発用）

MCP Inspectorを使うと、各サーバーのツールをWebインタフェースでテストできます：

```bash
# conda環境をアクティベート
conda activate mcp-md

# MCP Inspector起動（Structure Serverを例に）
mcp dev servers/structure_server.py

# 別のサーバーをテストする場合
mcp dev servers/genesis_server.py
mcp dev servers/complex_server.py
mcp dev servers/ligand_server.py
```

ブラウザが自動的に開き、以下が可能：
- 利用可能なツール一覧の表示
- 各ツールのスキーマ確認
- パラメータを入力してツールを実行
- レスポンスの確認

### ワークフロー例

#### 例1: Boltz-2でFASTAから構造予測 → MD入力

```python
# 1. 構造予測
result = await boltz2_protein_from_seq(
    sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL...",
    use_msa=True,
    num_models=5
)

# 2. 系の構築
system = await build_system_tleap(
    protein_pdb=result["structures"][0],
    forcefield="leaprc.protein.ff19SB"
)

# 3. エネルギー最小化
minimized = await openmm_minimize(
    prmtop=system["prmtop"],
    inpcrd=system["inpcrd"],
    max_iterations=5000
)
```

#### 例2: FASTA + SMILES → 複合体 + 親和性予測

```python
# タンパク質-リガンド複合体の構造と親和性を同時予測
result = await boltz2_complex(
    protein_fasta="MKTAYIAK...",
    ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    use_msa=True
)

# 親和性結果
print(f"Binder probability: {result['affinity']['probability_binary']:.2f}")
print(f"IC50: {result['affinity']['ic50_um']:.2f} μM")

# 複合体構造でMD系構築
complex_system = await build_system_tleap(
    protein_pdb=result["structures"][0],
    ligand_lib="ligand.lib",
    ligand_frcmod="ligand.frcmod"
)
```

#### 例3: PDB + smina docking → MD入力

```python
# 既存PDB構造にリガンドをドッキング
protein = await fetch_pdb("1ABC")
cleaned = await clean_structure(protein["file_path"])

# SMILESから3D構造生成 + GAFF2パラメータ化
ligand_params = await parameterize_ligand_complete(
    smiles="CC(=O)Oc1ccccc1C(=O)O",
    charge_method="bcc"  # AM1-BCC
)

# sminaでドッキング
docked = await smina_dock(
    receptor=cleaned["output"],
    ligand=ligand_params["gaff_mol2"],
    center_x=10.0, center_y=15.0, center_z=20.0,
    size_x=20.0, size_y=20.0, size_z=20.0
)

# 複合体MD系構築
final_system = await build_system_tleap(
    protein_pdb=cleaned["output"],
    ligand_lib=ligand_params["library"],
    ligand_frcmod=ligand_params["frcmod"]
)
```

## ディレクトリ構造

```
mcp-md/
├── servers/              # FastMCP サーバー実装（7サーバー）
│   ├── structure_server.py   # PDB取得・修復
│   ├── genesis_server.py     # Boltz-2構造生成
│   ├── complex_server.py     # 複合体予測・ドッキング
│   ├── ligand_server.py      # 配位子パラメータ化
│   ├── assembly_server.py    # 系の組立
│   ├── export_server.py      # 形式変換
│   └── qc_min_server.py      # 品質チェック・最小化
├── common/               # 共通ライブラリ
│   ├── base.py          # BaseToolWrapper（外部ツール実行）
│   └── utils.py         # 共通ユーティリティ関数
├── core/                 # LangGraphエージェント実装
│   ├── langgraph_agent.py    # LangGraph Agent + MCP Client
│   ├── workflow_graph.py     # StateGraph定義
│   ├── workflow_nodes.py     # ノード実装
│   ├── workflow_state.py     # WorkflowState定義
│   ├── mcp_integration.py    # langchain-mcp-adapters統合
│   ├── decision_logger.py    # 意思決定ログ
│   └── models.py             # Pydanticモデル
├── checkpoints/          # LangGraphチェックポイント
│   └── workflow.db      # SQLiteステート保存
├── tests/                # テストコード
├── .cursor/              # Cursorプロジェクト設定
│   └── rules/           # プロジェクトルール
│       ├── project-rules.md         # プロジェクト全体のルール
│       └── notebook-development.md  # Notebook開発ルール
├── pyproject.toml        # プロジェクト設定（langchain統合）
├── ARCHITECTURE.md       # 詳細アーキテクチャ・技術仕様
├── AGENTS.md             # Cursor AI Agent設定
└── README.md             # このファイル
```

## 開発

### テスト実行

```bash
conda activate mcp-md
pytest tests/
```

### コードフォーマット

```bash
conda activate mcp-md

# フォーマット適用
black servers/ core/ common/

# Lintチェック
ruff check servers/ core/ common/

# 型チェック
mypy servers/ core/ common/
```

## 開発ワークフロー

### 新しいMCPサーバーの追加（FastMCP）

1. **サーバーファイル作成** (`servers/`)

   ```python
   # servers/new_server.py
   from pathlib import Path
   from fastmcp import FastMCP
   from common.base import BaseToolWrapper
   from common.utils import setup_logger, ensure_directory
   
   logger = setup_logger(__name__)
   mcp = FastMCP("New Server")
   
   # 作業ディレクトリ
   WORKING_DIR = Path("output/new_server")
   ensure_directory(WORKING_DIR)
   
   # 外部ツールラッパー（必要に応じて）
   tool_wrapper = BaseToolWrapper("tool_name", conda_env="mcp-md")
   
   @mcp.tool
   def process_data(input_file: str, param: int = 0) -> dict:
       """Process data with new tool
       
       Args:
           input_file: Input file path
           param: Optional parameter
       
       Returns:
           Processing results
       """
       logger.info(f"Processing {input_file}")
       
       # 外部ツール実行
       result = tool_wrapper.run(['-i', input_file, '--param', str(param)])
       
       return {
           "status": "success",
           "output_file": str(WORKING_DIR / "output.dat")
       }
   
   if __name__ == "__main__":
       mcp.run()  # STDIO transport (default)
   ```

2. **テスト作成** (`tests/`)

   ```python
   # tests/test_new_server.py
   import pytest
   from fastmcp import Client
   
   @pytest.mark.asyncio
   async def test_new_server():
       # Import server module to get mcp instance
       from servers import new_server
       
       # Connect to server using in-memory transport
       async with Client(new_server.mcp) as client:
           tools = await client.list_tools()
           assert "process_data" in [t.name for t in tools]
           
           result = await client.call_tool("process_data", {
               "input_file": "test.dat"
           })
           assert result.content[0].text  # Check result exists
   ```

3. **LangGraphに登録** (`core/mcp_integration.py`)

   ```python
   # create_mcp_client() のserver_configに追加
   "new_server": {
       "transport": "stdio",
       "command": python_exe,
       "args": ["-m", "servers.new_server"]
   },
   ```

### MCPツールの追加（FastMCP）

既存サーバーに新しいツールを追加する場合：

1. `servers/xxx_server.py`に`@mcp.tool`デコレータで関数追加
2. 型ヒントとdocstringで自動的にスキーマ生成
3. テスト追加

**例**: Structure Serverに新しいツール追加

```python
# servers/structure_server.py

@mcp.tool
def analyze_structure(pdb_file: str, analysis_type: str = "basic") -> dict:
    """Perform structure analysis
    
    Args:
        pdb_file: Input PDB file path
        analysis_type: Type of analysis (basic, detailed, full)
    
    Returns:
        Analysis results with metrics
    """
    logger.info(f"Analyzing {pdb_file}: {analysis_type}")
    
    if not Path(pdb_file).is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    # 解析実装
    metrics = {
        "num_atoms": count_atoms_in_pdb(pdb_file),
        "chains": get_pdb_chains(pdb_file),
        "analysis_type": analysis_type
    }
    
    return {
        "status": "success",
        "metrics": metrics
    }
```

### デバッグ方法

#### MCPサーバーのデバッグ

```bash
# サーバーを直接実行（フォアグラウンド）
conda activate mcp-md
python -m servers.structure_server

# 詳細ログ有効化
export MCP_MD_LOG_LEVEL=DEBUG
python -m servers.structure_server
```

#### Pythonデバッガ使用

```python
# サーバーコード内にブレークポイント設定
import pdb; pdb.set_trace()

# または
breakpoint()  # Python 3.7+
```

### プロジェクト詳細情報

詳細なアーキテクチャ、FastMCP統合の実装状況、技術仕様は`ARCHITECTURE.md`を参照してください。

## サポートされる力場

- **タンパク質**: ff19SB, ff14SB
- **核酸**: OL15, OL3
- **脂質**: lipid17, CHARMM36
- **糖鎖**: GLYCAM06
- **小分子**: GAFF2, OpenFF (SMIRNOFF)

## サポートされるMDエンジン

- **OpenMM**: フルサポート（推奨）
- **Amber**: prmtop/inpcrd出力
- **GROMACS**: ParmEdで変換
- **CHARMM**: PSF/CRD変換
- **NAMD**: PSF/PDB + 設定ファイル

## ライセンス

MIT License

## 引用

このツールを使用する場合、以下を引用してください：

### Boltz-2

```
S. Passaro et al., Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction.
```

### AmberTools

```
D. A. Case et al., AmberTools, J. Chem. Inf. Model. 63, 6183 (2023).
```

### OpenMM

```
P. Eastman et al., OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials, J. Phys. Chem. B 128, 109 (2024).
```

## コントリビューション

Issue、Pull Requestを歓迎します。

## サポート

バグ報告や機能要望は、GitHubのIssueでお願いします。
