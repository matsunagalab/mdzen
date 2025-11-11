# MCP-MD: Molecular Dynamics Input File Generation Agent

Amber系に特化したMD入力ファイル生成AIエージェントシステム。LangGraph + FastMCPで構築された3フェーズワークフロー（Clarification → Setup → Validation）。

## 特徴

- **LangGraph統合**: ステートフルなワークフロー、永続化、人間フィードバック
  - LangChain 1.0準拠のStateGraphベースの実装
  - langchain-mcp-adaptersで公式MCP統合
  - チェックポイント機能で中断・再開可能
- **Boltz-2統合**: FASTAやSMILESから高精度な構造予測と結合親和性予測
- **AmberTools完結**: 配位子パラメータ化に外部QMソフト不要（AM1-BCC電荷計算）
- **FastMCP統合**: モジュラーな7つの独立サーバー、型安全な自動スキーマ生成
- **OpenMM専用**: Pythonプログラマブルなプロダクションレディなスクリプト生成

## 📚 ドキュメント

- **[ARCHITECTURE.md](ARCHITECTURE.md)** - プロジェクト全体のアーキテクチャ・実装プラン・技術仕様
- **[AGENTS.md](AGENTS.md)** - Cursor AI Agent設定とガイドライン
- **[.cursor/rules/](.cursor/rules/)** - プロジェクトルールと開発ワークフロー

## インストール

### 前提条件

- Python 3.11以上
- [conda](https://docs.conda.io/en/latest/) または [mamba](https://mamba.readthedocs.io/)
- GPU推奨（Boltz-2、OpenMM高速化）

### 手順

#### 1. conda環境のセットアップ

```bash
# conda環境作成
conda create -n mcp-md python=3.11
conda activate mcp-md

# 科学計算パッケージをインストール
conda install -c conda-forge openmm rdkit mdanalysis biopython pandas numpy scipy openblas pdbfixer

# MD準備ツール
conda install -c conda-forge ambertools packmol smina
```

#### 2. Pythonパッケージのインストール

```bash
# プロジェクトのクローン
git clone https://github.com/matsunagalab/mcp-md.git
cd mcp-md

# パッケージをインストール（editable mode）
pip install -e .
```

#### 3. Boltz-2のインストール（オプション）

Boltz-2は Phase 2-3（Setup/Validation）で使用します。必要になったときにインストールしてください：

```bash
# CUDA対応GPUがある場合
pip install 'boltz[cuda]' --no-deps

# その後、不足している依存関係を個別にインストール
pip install torch hydra-core pytorch-lightning einops einx mashumaro modelcif wandb

# または、scipyをダウングレードしてから通常インストール
conda install -c conda-forge scipy=1.13.1
pip install 'boltz[cuda]'
```

> **注意**: Boltz-2の依存関係の一つ（fairscale）がscipy==1.13.1を厳密に要求するため、condaで既にインストールされているscipyと競合する場合があります。`--no-deps`オプションを使用することで、既存のパッケージを保持したまま、不足しているものだけを追加できます。

#### 4. Ollamaのインストール（オプション）
OllamaはLocal LLMのローカル実行環境です。デフォルトではOllamaの`gpt-oss:20b`モデルを使用します。

```bash
# Macの場合
brew install ollama
brew pull gpt-oss:20b
brew services start ollama
```

## 使用方法

### Phase 1: Clarification（要件明確化とブリーフ生成）

```bash
# テストスクリプトを実行
python demo_clarification.py
```

開発は、Notebookで開発：

```bash
jupyter notebook
# notebooks/1_clarification.ipynb を開く
```

### MCPサーバーのテスト

各FastMCPサーバーを単独でテスト可能：

```bash
# MCP Inspector起動（Structure Serverを例に）
mcp dev servers/structure_server.py

# 別のサーバーをテストする場合
mcp dev servers/genesis_server.py
mcp dev servers/complex_server.py
mcp dev servers/ligand_server.py
```

## ディレクトリ構造

```
mcp-md/
├── notebooks/            # 🎯 開発の中心（Notebook-first development）
│   ├── 1_clarification.ipynb       # Phase 1: 要件明確化
│   ├── 2_setup_agent.ipynb         # Phase 2基本: Setupエージェント
│   ├── 3_setup_coordinator.ipynb   # Phase 2高度: Coordinator-Tools
│   ├── 4_validation.ipynb          # Phase 3: 検証・エクスポート
│   ├── 5_full_agent.ipynb          # 全統合: End-to-End
│   └── utils.py                    # Notebook用ユーティリティ
│
├── src/mcp_md/           # 生成されたソースコード（%%writefileで自動生成）
│   ├── state_scope.py              # Phase 1状態定義
│   ├── state_setup.py              # Phase 2状態定義
│   ├── clarification_agent.py      # Phase 1実装
│   ├── setup_agent.py              # Phase 2基本実装
│   ├── setup_coordinator.py        # Phase 2高度実装
│   ├── validation_agent.py         # Phase 3実装
│   ├── full_agent.py               # 統合実装
│   ├── prompts.py                  # プロンプトテンプレート
│   ├── mcp_integration.py          # MCP統合
│   └── utils.py                    # ユーティリティ
│
├── servers/              # FastMCPサーバー（7サーバー）
│   ├── structure_server.py         # PDB取得・修復
│   ├── genesis_server.py           # Boltz-2構造生成
│   ├── complex_server.py           # 複合体予測・ドッキング
│   ├── ligand_server.py            # 配位子パラメータ化
│   ├── assembly_server.py          # 系の組立
│   ├── export_server.py            # 形式変換
│   └── qc_min_server.py            # 品質チェック・最小化
│
├── common/               # 共通ライブラリ
│   ├── base.py                     # BaseToolWrapper
│   └── utils.py                    # 共通ユーティリティ
│
├── checkpoints/          # LangGraphチェックポイント
├── runs/                 # 実行結果
├── ARCHITECTURE.md       # 詳細アーキテクチャ
├── AGENTS.md             # Cursor AI Agent設定
└── README.md             # このファイル
```

## 開発ワークフロー

### Notebook-First Development

このプロジェクトは **Notebook-First Development** を採用しています：

```
✅ notebooks/*.ipynb を編集
✅ %%writefile で src/mcp_md/ を生成
✅ Notebookでテスト・実行

🚫 src/mcp_md/ を直接編集しない
```

詳細は [.cursor/rules/notebook-development.md](.cursor/rules/notebook-development.md) を参照。

### コードフォーマット

```bash
# フォーマットチェック
ruff check src/mcp_md/

# 自動修正
ruff check src/mcp_md/ --fix
```

> **注意**: フォーマット問題が見つかった場合、`src/`ファイルではなく、**Notebookの`%%writefile`セル**で修正してください。

## ライセンス

MIT License

## 引用

このツールを使用する場合、以下を引用してください：

### Boltz-2

```
S. Passaro et al., Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction.
bioRxiv (2025). doi:10.1101/2025.06.14.659707
```

### AmberTools

```
D. A. Case et al., AmberTools, J. Chem. Inf. Model. 63, 6183 (2023).
```

### OpenMM

```
P. Eastman et al., OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials,
J. Phys. Chem. B 128, 109 (2024).
```
