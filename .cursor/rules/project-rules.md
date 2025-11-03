# MCP-MD - プロジェクトルール

## プロジェクト概要

このリポジトリは、LangGraphを使用してAmber系MD準備に特化したAIエージェントシステムを構築します。
**deep_research_from_scratch**の実装パターンを完全に踏襲し、5つのチュートリアルNotebookを通じて段階的に実装します。

## リポジトリ構造

```
mcp-md/
├── notebooks/              # 🎯 開発の中心（MODIFY THESE）✏️
│   ├── 1_clarification.ipynb     # Phase 1: 要件明確化
│   ├── 2_setup_agent.ipynb       # Phase 2基本: Setupエージェント
│   ├── 3_setup_coordinator.ipynb # Phase 2高度: Coordinator-Tools
│   ├── 4_validation.ipynb        # Phase 3: 検証・エクスポート
│   ├── 5_full_agent.ipynb        # 全統合: End-to-End
│   └── utils.py                  # Notebook用ユーティリティ
│
├── src/mcp_md/            # 🚫 生成されたソースコード（DO NOT MODIFY）
│   ├── state_scope.py
│   ├── state_setup.py
│   ├── clarification_agent.py
│   ├── setup_agent.py
│   ├── setup_coordinator.py
│   ├── validation_agent.py
│   ├── full_agent.py
│   ├── prompts.py
│   ├── mcp_integration.py
│   ├── decision_logger.py
│   ├── report_generation.py
│   └── utils.py
│
├── servers/               # ✅ FastMCPサーバー（既存・維持）
│   ├── structure_server.py
│   ├── genesis_server.py
│   ├── complex_server.py
│   ├── ligand_server.py
│   ├── assembly_server.py
│   ├── export_server.py
│   └── qc_min_server.py
│
├── common/                # ✅ 共通ライブラリ（既存・維持）
│   ├── base.py
│   └── utils.py
│
├── ARCHITECTURE.md        # アーキテクチャドキュメント
└── AGENTS.md              # Cursor AI Agent設定
```

## システムアーキテクチャ

### 3フェーズワークフロー

1. **Clarification** (Notebook 1): 要件明確化と構造化ブリーフ生成
2. **Setup** (Notebooks 2-3): 固定スケルトンに沿ったMDシステム構築
3. **Validation & Export** (Notebook 4): QC検証とファイルエクスポート

### 主要コンポーネント

- **Clarification Agent**: ユーザー意図を明確化し、シミュレーションブリーフを生成
- **Setup Agent**: ReActパターンでMCPツールを使用してセットアップ実行
- **Setup Coordinator**: Supervisorパターンで高度なツール選択と決定ログ記録
- **Validation Agent**: QC検証、形式変換、レポート生成
- **Full System**: 全コンポーネントをエンドツーエンドワークフローに統合

## deep_research_from_scratchとの対応

### アーキテクチャの対応

| deep_research | mcp-md | 目的 |
|--------------|--------|------|
| **Phase 1: Scope** | **Phase 1: Clarification** | 要件明確化 |
| `clarify_with_user` | `clarify_requirements` | 明確化判定 |
| `write_research_brief` | `generate_simulation_brief` | ブリーフ生成 |
| **Phase 2: Research** | **Phase 2: Setup** | 実行フェーズ |
| `researcher` + `researcher_tools` | `setup_agent` + `setup_tools` | 基本エージェント |
| `supervisor` + `supervisor_tools` | `setup_coordinator` + `setup_tools` | Coordinatorパターン |
| **Phase 3: Write** | **Phase 3: Validation** | 結果生成 |

### Structured Outputの対応

| deep_research | mcp-md | 用途 |
|--------------|--------|------|
| `ClarifyWithUser` | `ClarifyWithUser` | 明確化判定 |
| `ResearchQuestion` | `SimulationBrief` | ブリーフ生成 |
| `ConductResearch` | `ExecuteSetupStep` | タスク委譲 |
| `ResearchComplete` | `SetupComplete` | 完了シグナル |

### 状態管理の対応

| deep_research | mcp-md | 用途 |
|--------------|--------|------|
| `AgentInputState` | `AgentInputState` | ユーザー入力 |
| `AgentState` | `AgentState` | メイン状態 |
| `SupervisorState` | `SetupState` | サブグラフ状態 |
| `SupervisorOutputState` | `SetupOutputState` | サブグラフ出力 |

## 重要な違い

### 並列 vs 順次実行

- **deep_research**: 並列研究（`asyncio.gather()`で複数のresearch agentを同時実行）
- **mcp-md**: **順次実行**（固定スケルトンの各ステップを順番に実行）

### ツール統合

- **deep_research**: Tavily検索 + MCP（ファイルシステム）
- **mcp-md**: FastMCPサーバー（7つの専門サーバー）

## 参考資料

### 必読
1. **deep_research_from_scratch**
   - README.md: 全体アーキテクチャ
   - CLAUDE.md: 開発ワークフロー
   - notebooks/: 5つの実装例

2. **LangGraph公式ドキュメント**
   - Command API
   - Subgraphs
   - Persistence

3. **MCP Protocol**
   - langchain-mcp-adapters
   - MCP Specification

