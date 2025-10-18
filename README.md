# MCP-MD: Molecular Dynamics Input File Generation Agent

CHARMM-GUIに代わる、お手軽でフレクシブルなMD入力ファイル生成システム。Boltz-2による構造・親和性予測、AmberToolsによる配位子パラメータ化、OpenMMによるMD実行を統合。

## 特徴

- **Boltz-2統合**: FASTAやSMILESから高精度な構造予測と結合親和性予測
  - 構造予測、複合体+親和性予測、バーチャルスクリーニング、欠損残基補完
- **AmberTools完結**: 配位子パラメータ化に外部QMソフト不要（AM1-BCC電荷計算）
  - SMILES → 3D → GAFF2パラメータ → tleapライブラリの完全自動化
- **高度な構造処理** 🆕:
  - PDB2PQR+PROPKAによるpH指定プロトネーション
  - ジスルフィド結合・金属サイト自動検出
- **膜タンパク質系** 🆕: Packmol-Memgen統合で脂質二重層自動構築
- **OpenMM専用**: Pythonプログラマブルなプロダクションレディなスクリプト生成
- **LM Studio統合**: ローカルLLMによる自然言語ワークフロー生成
- **MCPアーキテクチャ**: 機能別サーバーでモジュラーな設計

## 📚 ドキュメント

- **[Phase 1/2/4実装詳細](docs/PHASE_124_IMPLEMENTATION.md)** - 構造・配位子・系組立の完全実装ガイド
- **[Phase 1/2/4統合ワークフロー](examples/phase_124_workflow.md)** - 実践的な使用例とコード

## セットアップ

### 前提条件

- Python 3.11以上
- [uv](https://github.com/astral-sh/uv) (Pythonパッケージマネージャー)
- [conda](https://docs.conda.io/en/latest/) または [mamba](https://mamba.readthedocs.io/) (外部ツール用)
- [LM Studio](https://lmstudio.ai/) (ローカルLLM実行)
- GPU推奨（Boltz-2、OpenMM高速化）

### インストール手順

#### 1. リポジトリのクローン

```bash
git clone <repository-url>
cd mcp-md
```

#### 2. Python環境のセットアップ (uv)

```bash
# uv仮想環境作成
uv venv

# 仮想環境をアクティベート
source .venv/bin/activate  # Linux/macOS
# または
.venv\Scripts\activate  # Windows

# 依存関係をインストール
uv pip install -e .

# Boltz-2インストール（GPU版）
uv pip install "boltz[cuda]" -U

# 開発用パッケージ（オプション）
uv pip install -e ".[dev]"
```

#### 3. 外部ツールのインストール (conda)

```bash
# conda環境作成
conda create -n mcp-md-tools python=3.11
conda activate mcp-md-tools

# AmberTools, Packmol, smina
conda install -c conda-forge ambertools packmol smina

# PDB2PQR, PROPKA
pip install pdb2pqr propka
```

#### 4. LM Studioのセットアップ

1. [LM Studio](https://lmstudio.ai/)をダウンロード・インストール
2. LM Studio GUIでモデルをダウンロード（推奨: `gpt-oss-20b`）
3. `Local Server`タブで`Start Server`をクリック（デフォルト: `http://localhost:1234`）
4. 環境変数を設定:

```bash
export LM_STUDIO_BASE_URL="http://localhost:1234/v1"
export LM_STUDIO_MODEL="gpt-oss-20b"
```

## 使用方法

### MCPサーバーの起動

各機能は独立したMCPサーバーとして動作します：

```bash
# Structure Server（構造取得・Boltz-2予測）
uv run python servers/structure_server.py

# Ligand Server（配位子パラメータ化）
uv run python servers/ligand_server.py

# Docking Server（smina ドッキング）
uv run python servers/docking_server.py

# Assembly Server（系の組立）
uv run python servers/assembly_server.py

# Protocol Server（OpenMM MDスクリプト）
uv run python servers/protocol_server.py

# Export Server（形式変換）
uv run python servers/export_server.py
```

### ワークフロー例

#### 例1: Boltz-2でFASTAから構造予測 → MD入力

```python
# 1. 構造予測
result = await predict_structure_boltz2(
    fasta="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL...",
    use_msa=True,
    num_models=5
)

# 2. 系の構築
system = await build_protein_system(
    pdb_file=result["structures"][0],
    forcefield="ff19SB"
)

# 3. 溶媒化・イオン追加
solvated = await solvate_box(system, padding=10.0)
final = await add_ions(solvated, concentration=0.15)

# 4. OpenMM MDスクリプト生成
workflow = await create_openmm_workflow(
    prmtop=final["prmtop"],
    inpcrd=final["inpcrd"],
    protocol="standard"
)
```

#### 例2: FASTA + SMILES → 複合体 + 親和性予測

```python
# タンパク質-リガンド複合体の構造と親和性を同時予測
result = await predict_complex_with_affinity(
    protein_fasta="MKTAYIAK...",
    ligand_smiles=["CC(=O)Oc1ccccc1C(=O)O"],  # Aspirin
    predict_affinity=True
)

# 親和性結果
print(f"Binder probability: {result['affinity']['probability_binary']:.2f}")
print(f"IC50: {result['affinity']['ic50_um']:.2f} μM")

# 複合体構造でMD系構築
complex_system = await build_complex_system(
    protein_pdb=result["structures"][0],
    ligand_mol2="ligand.mol2"
)
```

#### 例3: PDB + smina docking → MD入力

```python
# 既存PDB構造にリガンドをドッキング
protein = await fetch_pdb("1ABC")
cleaned = await clean_structure(protein)

# SMILESから3D構造生成
ligand_3d = await smiles_to_3d("CC(=O)Oc1ccccc1C(=O)O")

# AmberToolsでGAFF2パラメータ生成
params = await generate_gaff_params(
    ligand_file=ligand_3d["mol2"],
    charge_method="bcc"  # AM1-BCC
)

# sminaでドッキング
docked = await dock_ligand_smina(
    receptor_pdb=cleaned["pdb"],
    ligand_mol2=params["mol2"],
    center=[10.0, 15.0, 20.0],
    size=[20.0, 20.0, 20.0]
)

# 複合体MD系構築
final_system = await build_complex_system(
    protein_pdb=cleaned["pdb"],
    ligand_mol2=docked["poses"][0]
)
```

## ディレクトリ構造

```
mcp-md/
├── servers/              # MCPサーバー実装
│   ├── structure_server.py
│   ├── ligand_server.py
│   ├── docking_server.py
│   ├── assembly_server.py
│   ├── protocol_server.py
│   └── export_server.py
├── core/                 # コアエンジン
│   ├── llm_client.py     # LM Studio統合
│   ├── planner.py        # ワークフロープランニング
│   ├── validator.py      # QC/バリデーション
│   └── workflow.py       # ワークフロー実行
├── tools/                # 外部ツールラッパー
│   ├── boltz2_wrapper.py
│   ├── pdbfixer_wrapper.py
│   ├── openmm_wrapper.py
│   ├── rdkit_wrapper.py
│   ├── ambertools_wrapper.py
│   ├── smina_wrapper.py
│   └── packmol_wrapper.py
├── tests/                # テストコード
├── examples/             # 使用例・YAMLテンプレート
├── docs/                 # ドキュメント
├── pyproject.toml        # プロジェクト設定
└── README.md             # このファイル
```

## 開発

### テスト実行

```bash
uv run pytest tests/
```

### コードフォーマット

```bash
# フォーマット適用
uv run black servers/ core/ tools/

# Lintチェック
uv run ruff check servers/ core/ tools/

# 型チェック
uv run mypy servers/ core/ tools/
```

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

```bibtex
@article{passaro2025boltz2,
  author = {Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and ...},
  title = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
  year = {2025},
  journal = {bioRxiv}
}
```

### AmberTools

```bibtex
@article{case2023ambertools,
  title={AmberTools},
  author={Case, D.A. and ...},
  journal={Journal of Chemical Information and Modeling},
  year={2023}
}
```

### OpenMM

```bibtex
@article{eastman2017openmm,
  title={OpenMM 7: Rapid development of high performance algorithms for molecular dynamics},
  author={Eastman, Peter and ...},
  journal={PLOS Computational Biology},
  year={2017}
}
```

## コントリビューション

Issue、Pull Requestを歓迎します。

## サポート

バグ報告や機能要望は、GitHubのIssueでお願いします。

