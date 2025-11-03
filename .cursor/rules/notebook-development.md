# Notebook開発ルール

## 🚨 最重要ルール

**Notebooksが真実の源泉（Source of Truth）です。このルールは絶対です。**

```
notebooks/     ← MODIFY THESE（開発の中心）✏️
src/mcp_md/    ← GENERATED CODE（直接編集禁止）🚫
```

`src/mcp_md/`のソースコードは、Notebooksから`%%writefile`マジックコマンドを使って自動生成されます。

## コード生成の仕組み

### 1. Notebooksに`%%writefile`セルが含まれる

各NotebookはJupyterの`%%writefile`マジックを使って、`src/`にコードを直接書き出します。

```python
%%writefile ../src/mcp_md/clarification_agent.py

"""
Clarification Agent Implementation
"""
from langchain.chat_models import init_chat_model
from langgraph.graph import StateGraph, START, END
# ... 実装コード ...
```

### 2. Notebooksは実行可能なチュートリアル

概念をインタラクティブに示しながら、本番コードを生成します。

```python
# Notebookセル1: コード生成
%%writefile ../src/mcp_md/clarification_agent.py
# ... 実装コード ...

# Notebookセル2: 即座にテスト
from mcp_md.clarification_agent import clarify_requirements

result = await clarify_requirements(test_state)
print(result)
```

### 3. ソースファイルは生成物

`src/`の`.py`ファイルは出力であり、入力ではありません。

## 開発ガイドライン

### ✅ DO（推奨）

- **DO**: `notebooks/`ディレクトリのNotebooksを編集
- **DO**: Notebookセルを実行してソースコードを再生成
- **DO**: Notebooksを実行して変更をテスト
- **DO**: `%%writefile`で生成されたコードをその場でテスト
- **DO**: Notebookに説明とコードを一緒に記述

### ❌ DON'T（禁止）

- **DON'T**: `src/mcp_md/`のファイルを直接編集
- **DON'T**: `src/`への手動変更が永続化されることを期待
- **DON'T**: `src/`ファイルを手動で作成
- **DON'T**: Notebookを経由せずにソースコードを修正

## 開発サイクル

### 1. Notebookで実装

```python
%%writefile ../src/mcp_md/clarification_agent.py

"""Clarification Agent Implementation"""

from langchain.chat_models import init_chat_model
# ... 実装コード ...
```

### 2. その場でテスト（同じNotebook内）

```python
# 次のセルで即座にテスト
from mcp_md.clarification_agent import clarify_requirements

result = await clarify_requirements(test_state)
print(result)
```

### 3. 修正が必要な場合

- ❌ `src/mcp_md/clarification_agent.py`を直接編集 → ダメ！
- ✅ Notebookの`%%writefile`セルに戻る → 修正 → 再実行

## 実装の進め方

### Notebook 1: Clarification (3-4日)
- Structured Output定義（ClarifyWithUser, SimulationBrief）
- Command APIでルーティング
- サブグラフ構築

**生成ファイル**:
- `src/mcp_md/state_scope.py`
- `src/mcp_md/clarification_agent.py`
- `src/mcp_md/prompts.py`（初版）

### Notebook 2: Setup Agent (3-4日)
- MCP統合（MultiServerMCPClient）
- 固定スケルトン実装
- ReActエージェント

**生成ファイル**:
- `src/mcp_md/state_setup.py`
- `src/mcp_md/setup_agent.py`
- `src/mcp_md/mcp_integration.py`
- `src/mcp_md/utils.py`

### Notebook 3: Setup Coordinator (4-5日)
- Coordinator-Toolsパターン
- Structured Tools（ExecuteSetupStep, SetupComplete）
- 決定ログ記録

**生成ファイル**:
- `src/mcp_md/setup_coordinator.py`
- `src/mcp_md/prompts.py`（Coordinatorプロンプト追加）
- `src/mcp_md/decision_logger.py`

### Notebook 4: Validation (2-3日)
- QC検証
- 形式変換
- レポート生成

**生成ファイル**:
- `src/mcp_md/validation_agent.py`
- `src/mcp_md/report_generation.py`
- `src/mcp_md/prompts.py`（レポートプロンプト追加）

### Notebook 5: Full Agent (3-4日)
- 3フェーズ統合
- チェックポイント機能
- エンドツーエンドテスト

**生成ファイル**:
- `src/mcp_md/full_agent.py`
- `src/mcp_md/__init__.py`

## コード品質とフォーマット

### Ruffフォーマットチェック

生成されたソースファイル全体で一貫したコードフォーマットを維持するため、定期的にruffを実行：

```bash
# フォーマット問題をチェック
ruff check src/

# 可能な場合、フォーマット問題を自動修正
ruff check src/ --fix

# 特定ファイルをチェック
ruff check src/mcp_md/clarification_agent.py
```

**重要**: `src/`のソースファイルはNotebooksから生成されるため、フォーマット問題はソースファイルではなく、**Notebookの`%%writefile`セル**で修正する必要があります。Notebooksでフォーマットを修正した後、Notebookセルを実行してソースファイルを再生成してください。

### よくあるフォーマット修正

- **D212**: docstring要約が三重引用符と同じ行から始まることを確認
- **I001**: インポートを適切に整理（標準ライブラリ → サードパーティ → ローカルインポート）
- **F401**: 未使用のインポートを削除
- **D415**: docstring要約にピリオドを追加

## メリット

- **インタラクティブな開発**: Notebookで即座にテスト
- **バージョン管理がシンプル**: Notebookが真実の源泉
- **ドキュメント統合**: 説明とコードが一体
- **再現性**: Notebookを実行すればsrc/が生成される

