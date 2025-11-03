# MCP-MD - Cursor AI Agent設定

このファイルは、Cursor AI Agentがmcp-mdプロジェクトで作業する際のガイドラインを定義します。

## エージェントの役割

Cursor AI Agentは、deep_research_from_scratchの実装パターンに従って、LangGraphベースのMD準備エージェントシステムを構築する開発アシスタントです。

## 重要な原則

### 1. Notebookファースト開発（最重要）

**絶対ルール**: `notebooks/`のファイルのみを編集し、`src/mcp_md/`は直接編集しない。

```
✅ notebooks/*.ipynb を編集
✅ %%writefile で src/mcp_md/ を生成
✅ Notebookでテスト・実行

🚫 src/mcp_md/ を直接編集
🚫 手動で src/ ファイルを作成
```

### 2. deep_research_from_scratchパターン準拠

以下のパターンを完全に踏襲：

- **状態管理**: InputState / MainState / SubgraphState / OutputState
- **Structured Output**: Pydanticスキーマで決定を明示化
- **Command API**: ノード内での条件分岐とルーティング
- **Supervisor パターン**: coordinator + tools の2ノード構成

### 3. 段階的実装

5つのNotebookを順番に実装：

1. `1_clarification.ipynb` - Phase 1: Scoping
2. `2_setup_agent.ipynb` - Phase 2: Setup (基本)
3. `3_setup_coordinator.ipynb` - Phase 2: Setup (高度)
4. `4_validation.ipynb` - Phase 3: Validation
5. `5_full_agent.ipynb` - 全統合

## ファイル編集のルール

### 編集可能なファイル ✅

- `notebooks/*.ipynb` - メインの開発対象
- `notebooks/utils.py` - Notebook用ユーティリティ
- `servers/*.py` - FastMCPサーバー（既存維持）
- `common/*.py` - 共通ライブラリ（既存維持）
- `ARCHITECTURE.md` - アーキテクチャドキュメント
- `README.md` - プロジェクト説明

### 編集禁止のファイル 🚫

- `src/mcp_md/*.py` - Notebookから自動生成される
- `deep_research_from_scratch/` - 参考実装（読み取り専用）

## コード生成パターン

### Notebookでの実装パターン

```python
# ===== セル1: コード生成 =====
%%writefile ../src/mcp_md/clarification_agent.py

"""Clarification Agent Implementation"""

from langchain.chat_models import init_chat_model
from langgraph.graph import StateGraph, START, END
from langgraph.types import Command
# ... 実装 ...

# ===== セル2: 即座にテスト =====
from mcp_md.clarification_agent import clarify_requirements

# テストコード
result = await clarify_requirements(test_state)
print(result)
```

### %%writefileの使い方

- ファイルパスは `../src/mcp_md/` からの相対パス
- 必ずdocstringを含める（D212形式）
- インポートは標準ライブラリ → サードパーティ → ローカルの順

## Structured Outputパターン

### 明確化判定スキーマ

```python
class ClarifyWithUser(BaseModel):
    """ユーザー明確化の決定スキーマ"""
    need_clarification: bool = Field(description="追加の明確化が必要かどうか")
    question: str = Field(description="ユーザーに尋ねる具体的な質問")
    verification: str = Field(description="情報収集完了後の確認メッセージ")
```

### ブリーフ生成スキーマ

```python
class SimulationBrief(BaseModel):
    """シミュレーション要件の構造化スキーマ"""
    pdb_id: Optional[str] = Field(default=None, description="PDB ID")
    fasta_sequence: Optional[str] = Field(default=None, description="FASTA配列")
    ligand_smiles: str = Field(description="リガンドのSMILES")
    # ... その他パラメータ
```

## Command APIパターン

### ノード内ルーティング

```python
def clarify_requirements(state: AgentState) -> Command[Literal["generate_simulation_brief", "__end__"]]:
    """ユーザー要件の明確化ノード"""
    model_with_structured_output = model.with_structured_output(ClarifyWithUser)
    response = model_with_structured_output.invoke([...])
    
    if response.need_clarification:
        return Command(goto=END, update={"messages": [AIMessage(content=response.question)]})
    else:
        return Command(goto="generate_simulation_brief", update={"messages": [...]})
```

## Coordinator-Toolsパターン

### 2ノード構成

```python
# Coordinator ノード（決定）
async def setup_coordinator(state: SetupState) -> Command[Literal["setup_tools"]]:
    setup_tools = [ExecuteSetupStep, SetupComplete, think_tool]
    model_with_tools = model.bind_tools(setup_tools)
    response = await model_with_tools.ainvoke([...])
    return Command(goto="setup_tools", update={"setup_messages": [response]})

# Tools ノード（実行）
async def setup_tools(state: SetupState) -> Command[Literal["setup_coordinator", "__end__"]]:
    # ツール実行と決定ログ記録
    # ...
    return Command(goto="setup_coordinator", update={...})
```

## 品質保証

### Ruffチェック

コード生成後、必ずruffでフォーマットチェック：

```bash
ruff check src/mcp_md/
```

フォーマット問題が見つかった場合、**Notebookの`%%writefile`セルで修正**し、再実行。

### テスト

各Notebookで実装後、必ず動作確認：

1. Notebookセルを順番に実行
2. 生成されたコードをインポート
3. テスト入力で動作確認
4. 出力を検証

## エラーハンドリング

### よくあるエラーと対処法

#### 1. ImportError
```
ModuleNotFoundError: No module named 'mcp_md'
```

**対処**: `%%writefile`セルを再実行してファイルを生成

#### 2. フォーマットエラー
```
ruff: D212 Multi-line docstring summary should start at the first line
```

**対処**: Notebookの`%%writefile`セルでdocstringを修正し、再実行

#### 3. 状態管理エラー
```
KeyError: 'simulation_brief'
```

**対処**: 状態定義を確認し、必要なフィールドが初期化されているか確認

## 参考実装の活用

### deep_research_from_scratchを参照

プロジェクト内の`deep_research_from_scratch/`ディレクトリに参考実装があります：

- `notebooks/1_scoping.ipynb` → `1_clarification.ipynb`の参考
- `notebooks/4_research_supervisor.ipynb` → `3_setup_coordinator.ipynb`の参考
- `src/deep_research_from_scratch/` → 生成コードの参考

コピーペーストではなく、パターンを理解して適用してください。

## 開発の優先順位

1. **最優先**: Notebook 1（Clarification）の実装
2. **次**: Notebook 2（Setup Agent基本）の実装
3. **その後**: Notebook 3-5の順次実装
4. **最後**: 統合テストとドキュメント整備

## コミュニケーション

### ユーザーへの報告

- 各Notebookの完成時に進捗を報告
- 問題が発生した場合、deep_researchパターンと照らし合わせて説明
- 新しいファイルを生成した際は、`%%writefile`で生成されたことを明記

### 質問すべき事項

- deep_researchパターンからの逸脱が必要な場合
- MCPサーバーの動作に関する不明点
- テスト結果の解釈が必要な場合

## まとめ

**Cursor AI Agentとして最も重要なこと**:

1. ✅ **Notebooksを編集**し、`%%writefile`でソースを生成
2. ✅ **deep_researchパターン**を完全に踏襲
3. ✅ **段階的に実装**し、各Notebookで動作確認
4. 🚫 **src/を直接編集しない**

