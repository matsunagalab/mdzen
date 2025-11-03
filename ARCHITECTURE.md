# MCP-MD アーキテクチャ・実装プラン

## 1. プロジェクト概要

### 目的とポジショニング

**Amber系に最適化したAIエージェント＋MCPツール群**

- **主軸**: Amber/GAFF/OpenFF/ParmEd/OpenMM エコシステムに特化
- **非競合**: CHARMM-GUIとは棲み分け（CHARMM系は変換経由で二次対応、将来拡張）
- **永続化**: MCP標準でツール接続を維持可能（将来のLLM/実行基盤の更新に強い）
- **ホスト/クライアント**: [LangChain](https://github.com/langchain-ai/langchain) + [LangGraph](https://github.com/langchain-ai/langgraph)に統一（MCPツール統合）
- **参考実装**: [deep_research_from_scratch](https://github.com/langchain-ai/deep_research_from_scratch) の3フェーズアーキテクチャを適用

### 主要技術スタック

- **LangChain 1.0+**: LLM統合、ツール抽象化、プロンプト管理
  - LangChain 1.0では全てのchainsとagentsがLangGraph上に統一
  - `langchain-core`, `langchain-openai` (or `langchain-anthropic`)を使用
- **LangGraph 1.0+**: ステートフルなグラフベースのワークフロー構築
  - **Command API**: `Command(goto=..., update=...)`でノード内条件分岐
  - **Structured Output**: Pydanticモデルで決定を明示化・決定論化
  - **サブグラフ**: 各フェーズを独立したサブグラフとして実装
  - チェックポイント機能で永続化とtime-travel可能
  - 複雑な制御フローと人間フィードバック統合をネイティブサポート
- **FastMCP**: MCPサーバーの実装とツール提供
- **Boltz-2**: 構造予測・複合体生成ツール
- **AmberTools**: 完全OSS、配位子パラメータ化（GAFF2 + AM1-BCC）
- **OpenMM**: Pythonプログラマブル、GPU最適化、プロダクション対応MD
- **MCP (Model Context Protocol)**: 標準化されたツール統合（ツールの永続性・相互運用性）

### 主要機能

1. **3フェーズワークフロー**: Clarification → Setup → Validation & Export
2. **構造化意思決定**: Pydantic Structured Outputで決定を明示化
3. **品質保証**: 自作MolProbity等による物理化学的一貫性チェック
4. **再現性**: Plan/決定/生成物をJSON保存

---

## 2. 全体アーキテクチャ（3フェーズ設計）

### 3フェーズワークフロー（deep_research_from_scratchパターン）

```
┌─────────────────────────────────────────────────────────────────┐
│                   MCP-MD Agent System                            │
│                                                                   │
│  Phase 1: Clarification (Scope)                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ clarify_requirements → generate_simulation_brief          │   │
│  │                                                            │   │
│  │ Input:  User query ("PDBにAspirinをドッキング")           │   │
│  │ Output: SimulationBrief (structured)                      │   │
│  │         - pdb_id or fasta_sequence                        │   │
│  │         - ligand_smiles                                   │   │
│  │         - simulation_params (pH, salt, box, etc.)        │   │
│  │         - workflow_preferences                            │   │
│  └──────────────────────────────────────────────────────────┘   │
│                            ↓                                      │
│  Phase 2: Setup (Execute)                                        │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ setup_coordinator → setup_tools → [QC check]             │   │
│  │                                                            │   │
│  │ Subgraph with fixed skeleton:                            │   │
│  │   1. structure_fetch (PDB/Boltz-2)                       │   │
│  │   2. structure_repair (PDBFixer/PDB2PQR)                 │   │
│  │   3. ligand_param (GAFF2/AM1-BCC)                        │   │
│  │   4. complex_generation (Boltz-2/Smina)                  │   │
│  │   5. assembly (tleap, solvate, ions)                     │   │
│  │   6. qc_check (clash, bond, minimize)                    │   │
│  │                                                            │   │
│  │ Each step: Tool selection + Execution + Decision logging │   │
│  └──────────────────────────────────────────────────────────┘   │
│                            ↓                                      │
│  Phase 3: Validation & Export                                    │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ validate_system → export_files → generate_report         │   │
│  │                                                            │   │
│  │ Output: Final package                                     │   │
│  │   - prmtop, inpcrd (Amber)                               │   │
│  │   - Optional: GROMACS, OpenMM formats                    │   │
│  │   - qc_report.json                                       │   │
│  │   - decision_log.json                                    │   │
│  │   - metadata.json (再現性)                               │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘

[FastMCP Servers] (7 servers)
  ├─ Structure Server   (fetch, clean, protonate)
  ├─ Genesis Server     (Boltz-2 protein generation)
  ├─ Complex Server     (Boltz-2 complex, Smina dock)
  ├─ Ligand Server      (GAFF2/AM1-BCC parameterization)
  ├─ Assembly Server    (tleap, membrane, solvation)
  ├─ Export Server      (format conversion, packaging)
  └─ QC/Min Server      (minimization, validation)

[Persistent Storage]
  ├─ checkpoints/       (LangGraph state snapshots)
  │   └─ <thread_id>/   (会話セッション単位)
  └─ runs/<timestamp>/
      ├─ simulation_brief.json  (Phase 1 output)
      ├─ decision_log.json      (Phase 2 decisions)
      ├─ outputs/               (PDB, prmtop, inpcrd, etc.)
      ├─ qc_report.json         (Phase 3 validation)
      └─ metadata.json          (再現性情報)
```

### FastMCP統合の特徴

- **モジュラー設計**: 各サーバーファイルが完全に独立して動作可能
- **自動スキーマ生成**: 型ヒントとdocstringから自動的にMCPツールスキーマを生成
- **標準準拠**: MCP標準プロトコルに完全準拠
- **開発効率**: デコレータベースのシンプルなAPI（`@mcp.tool`）
- **独立実行**: 各サーバーが `python -m servers.{server_name}` で単独起動可能
- **共通ライブラリ**: `common/` モジュールで外部ツール実行とユーティリティを共有
- **LangChain統合**: `langchain-mcp-adapters`でMCPツールをLangChainツールとして利用

### 3フェーズ設計の詳細

#### Phase 1: Clarification (ユーザー要件の明確化)

**目的**: ユーザーの曖昧な要求を構造化されたシミュレーションブリーフに変換

**ノード構成**:
```python
clarify_requirements → generate_simulation_brief
```

**Structured Output**（Pydantic Schema）:
```python
class ClarifyWithUser(BaseModel):
    """ユーザー明確化の決定スキーマ"""
    need_clarification: bool = Field(
        description="追加の明確化が必要かどうか"
    )
    question: str = Field(
        description="ユーザーに尋ねる具体的な質問"
    )
    verification: str = Field(
        description="情報収集完了後の確認メッセージ"
    )

class SimulationBrief(BaseModel):
    """シミュレーション要件の構造化スキーマ"""
    pdb_id: Optional[str] = Field(description="PDB ID（既存構造の場合）")
    fasta_sequence: Optional[str] = Field(description="FASTA配列（de novo生成の場合）")
    ligand_smiles: str = Field(description="リガンドのSMILES文字列")
    
    # シミュレーションパラメータ
    ph: float = Field(default=7.4, description="pH値")
    salt_concentration: float = Field(default=0.15, description="塩濃度 (M)")
    water_model: str = Field(default="TIP3P", description="水モデル")
    box_padding: float = Field(default=12.0, description="Box padding (Å)")
    force_field: str = Field(default="ff19SB", description="タンパク質力場")
    
    # ワークフロー設定
    use_boltz2_docking: bool = Field(default=True, description="Boltz-2でドッキング")
    refine_with_smina: bool = Field(default=False, description="Sminaで精密化")
    output_formats: list[str] = Field(default=["amber"], description="出力形式")
```

**実装例**:
```python
def clarify_requirements(state: AgentState) -> Command[Literal["generate_simulation_brief", "__end__"]]:
    """ユーザー要件の明確化ノード"""
    model_with_structured_output = model.with_structured_output(ClarifyWithUser)
    
    response = model_with_structured_output.invoke([
        HumanMessage(content=clarify_prompt.format(
            messages=get_buffer_string(state["messages"]),
            date=get_today_str()
        ))
    ])
    
    if response.need_clarification:
        # 追加の質問が必要
        return Command(
            goto=END,
            update={"messages": [AIMessage(content=response.question)]}
        )
    else:
        # 情報十分、次のステップへ
        return Command(
            goto="generate_simulation_brief",
            update={"messages": [AIMessage(content=response.verification)]}
        )

def generate_simulation_brief(state: AgentState):
    """構造化されたシミュレーションブリーフを生成"""
    model_with_structured_output = model.with_structured_output(SimulationBrief)
    
    brief = model_with_structured_output.invoke([
        HumanMessage(content=brief_generation_prompt.format(
            messages=get_buffer_string(state["messages"]),
            date=get_today_str()
        ))
    ])
    
    return {
        "simulation_brief": brief,
        "setup_messages": [HumanMessage(content=f"Starting setup with: {brief.model_dump_json()}")]
    }
```

#### Phase 2: Setup (セットアップ実行)

**目的**: 固定スケルトンに沿ってMDシステムを構築、各ステップで最適なツールを選択

**Coordinator-Tools パターン**（deep_research supervisorパターン適用）:
```python
setup_coordinator → setup_tools → [next step or retry]
```

**固定スケルトン** (順序保証):
1. `structure_fetch` - PDB取得 or Boltz-2生成
2. `structure_repair` - PDBFixer + PDB2PQR
3. `ligand_param` - GAFF2/AM1-BCC
4. `complex_generation` - Boltz-2 or Smina
5. `assembly` - tleap系構築
6. `qc_check` - 品質チェック

**Structured Toolsで意思決定**:
```python
@tool
class ExecuteSetupStep(BaseModel):
    """セットアップステップ実行ツール"""
    step_name: str = Field(description="実行するステップ名")
    tool_name: str = Field(description="使用するツール名")
    parameters: dict = Field(description="ツールパラメータ")
    reason: str = Field(description="このツール選択の理由")

@tool
class SetupComplete(BaseModel):
    """セットアップ完了を示すツール"""
    pass

async def setup_coordinator(state: SetupState) -> Command[Literal["setup_tools"]]:
    """セットアップステップのコーディネータ"""
    setup_tools = [ExecuteSetupStep, SetupComplete, think_tool]
    model_with_tools = model.bind_tools(setup_tools)
    
    system_prompt = setup_coordinator_prompt.format(
        current_step=state["current_step"],
        simulation_brief=state["simulation_brief"],
        available_tools=get_available_tools_for_step(state["current_step"])
    )
    
    response = await model_with_tools.ainvoke(
        [SystemMessage(content=system_prompt)] + state["setup_messages"]
    )
    
    return Command(
        goto="setup_tools",
        update={"setup_messages": [response]}
    )

async def setup_tools(state: SetupState) -> Command[Literal["setup_coordinator", "__end__"]]:
    """セットアップツールの実行"""
    most_recent_message = state["setup_messages"][-1]
    
    if not most_recent_message.tool_calls:
        return Command(goto="setup_coordinator")
    
    tool_results = []
    decision_logs = []
    
    for tool_call in most_recent_message.tool_calls:
        if tool_call["name"] == "ExecuteSetupStep":
            # MCPツールを実行
            result = await execute_mcp_tool(
                tool_call["args"]["tool_name"],
                tool_call["args"]["parameters"]
            )
            tool_results.append(ToolMessage(content=result, tool_call_id=tool_call["id"]))
            
            # 決定ログに記録
            decision_logs.append({
                "step": tool_call["args"]["step_name"],
                "tool": tool_call["args"]["tool_name"],
                "parameters": tool_call["args"]["parameters"],
                "reason": tool_call["args"]["reason"],
                "timestamp": datetime.now().isoformat()
            })
        
        elif tool_call["name"] == "SetupComplete":
            # セットアップ完了
            return Command(
                goto=END,
                update={
                    "setup_messages": tool_results,
                    "decision_log": decision_logs
                }
            )
    
    return Command(
        goto="setup_coordinator",
        update={
            "setup_messages": tool_results,
            "decision_log": decision_logs
        }
    )
```

**特徴**:
- 固定スケルトンで順序保証（再現性）
- 各ステップで最適なツールを動的選択（柔軟性）
- すべての決定をStructured Outputで明示化
- 決定理由を必須化（説明可能性）

#### Phase 3: Validation & Export (検証とエクスポート)

**目的**: QC検証、形式変換、最終レポート生成

**ノード構成**:
```python
validate_system → export_files → generate_report
```

**実装例**:
```python
async def validate_system(state: AgentState):
    """システム全体の検証"""
    qc_results = await run_full_qc(
        prmtop=state["outputs"]["prmtop"],
        inpcrd=state["outputs"]["inpcrd"]
    )
    
    return {
        "qc_results": qc_results,
        "validation_passed": qc_results["overall_status"] == "pass"
    }

async def export_files(state: AgentState):
    """指定された形式でファイルをエクスポート"""
    exports = {}
    
    for format in state["simulation_brief"]["output_formats"]:
        if format == "amber":
            exports["amber"] = state["outputs"]
        elif format == "gromacs":
            exports["gromacs"] = await convert_to_gromacs(state["outputs"])
        elif format == "openmm":
            exports["openmm"] = await convert_to_openmm(state["outputs"])
    
    # パッケージング
    package_path = await package_system(exports, state["decision_log"], state["qc_results"])
    
    return {"exports": exports, "package_path": package_path}

async def generate_report(state: AgentState):
    """最終レポート生成"""
    report = {
        "simulation_brief": state["simulation_brief"].model_dump(),
        "decision_log": state["decision_log"],
        "qc_results": state["qc_results"],
        "outputs": state["exports"],
        "package_path": state["package_path"],
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "mcp_md_version": "0.1.0"
        }
    }
    
    # JSON保存
    report_path = save_report(report)
    
    # 人間可読なMarkdownレポート生成
    markdown_report = await generate_markdown_report(report)
    
    return {
        "final_report": markdown_report,
        "report_path": report_path,
        "messages": [AIMessage(content=f"Setup complete! Report: {report_path}")]
    }
```

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

## 4. LangGraph 1.0実装パターン

### LangChain 1.0とLangGraph 1.0の関係

- **LangChain 1.0の変更**: 従来の`chains`と`agents`を廃止、全てLangGraph上に統一
- **推奨アプローチ**: 
  - シンプルなReActエージェント → `create_react_agent()` (高レベル抽象化)
  - 複雑なワークフロー → LangGraphのStateGraphを直接使用（推奨）
- **本プロジェクトの選択**: 3フェーズワークフロー + 固定スケルトンのため、StateGraphを直接使用

### LangGraph 1.0の新機能

- **公式**: https://github.com/langchain-ai/langgraph
- **LangGraph 1.0の主要機能**: 
  - **Command API**: ノード内での条件分岐とルーティング（`Command(goto=..., update=...)`）
  - **Structured Output統合**: `model.with_structured_output()`で決定を明示化
  - **サブグラフ**: 各フェーズを独立したサブグラフとして実装、メイングラフで統合
  - **チェックポイント機能**: 永続化、time-travel、分岐実行
  - **Interrupt機能**: `interrupt_before/after`で人間フィードバック統合
- **LangChain統合**: LangChain `Tool`をノード内で直接利用可能
- **MCP統合**: `langchain-mcp-adapters`でMCPサーバーをLangChain ツールとして統合

### deep_research_from_scratchパターンの適用

我々のMCP-MDプロジェクトは、deep_research_from_scratchの以下のパターンを適用します：

| Deep Research | MCP-MD | 説明 |
|--------------|--------|------|
| **Scope Phase** | **Clarification Phase** | ユーザー要求の明確化と構造化ブリーフ生成 |
| `clarify_with_user` | `clarify_requirements` | Structured Outputで明確化の要否判定 |
| `write_research_brief` | `generate_simulation_brief` | 会話履歴から構造化ブリーフを生成 |
| **Research Phase** | **Setup Phase** | 研究実行 / MDセットアップ実行 |
| `supervisor` + `supervisor_tools` | `setup_coordinator` + `setup_tools` | Coordinator-Toolsパターン |
| Parallel research agents | Sequential setup steps | 並列研究 / 固定スケルトン |
| **Write Phase** | **Validation & Export Phase** | 結果統合 / 検証とエクスポート |
| `final_report_generation` | `validate_system` + `generate_report` | 最終レポート生成 |

### 運用のキーポイント

#### 1. 状態定義（deep_researchパターン完全適用）

**deep_researchの状態管理パターン**:
- `InputState`: ユーザー入力のみ（MessagesStateを継承）
- `MainState`: 全フィールド（MessagesStateを継承）
- `SubgraphState`: サブグラフ専用状態（TypedDict）
- `SubgraphOutputState`: サブグラフ出力状態（親へ返却）

**mcp-mdへの適用**:

```python
# src/mcp_md/state_scope.py (Notebook 1で生成)
import operator
from typing import TypedDict, Annotated, Sequence, Optional
from langgraph.graph import MessagesState, StateGraph, START, END
from langgraph.graph.message import add_messages
from langchain_core.messages import BaseMessage
from pydantic import BaseModel, Field

# ===== 入力状態（ユーザー入力のみ）=====
class AgentInputState(MessagesState):
    """
    エージェント入力状態 - ユーザーからのメッセージのみ.
    
    deep_researchの AgentInputState に対応。
    メイングラフの input_schema として使用。
    """
    pass

# ===== メイン状態（全フィールド）=====
class AgentState(MessagesState):
    """
    メインエージェント状態.
    
    deep_researchの AgentState に対応。
    MessagesStateを継承し、3フェーズ全体で使用されるフィールドを追加。
    
    注意: 一部のフィールドはサブグラフ状態と重複するが、これは
    deep_researchパターンに従った設計（状態の明示的な受け渡し）。
    """
    # Phase 1: Clarification
    research_brief: Optional[str] = None  # deep_researchと同じフィールド名
    simulation_brief: Optional[SimulationBrief] = None  # 構造化版
    
    # Phase 2: Setup
    setup_messages: Annotated[Sequence[BaseMessage], add_messages] = []
    decision_log: Annotated[list[dict], operator.add] = []
    outputs: dict = {}
    current_step: str = "structure_fetch"
    
    # Phase 3: Validation & Export
    qc_results: dict = {}
    exports: dict = {}
    package_path: str = ""
    final_report: str = ""

# ===== Structured Output スキーマ =====
class ClarifyWithUser(BaseModel):
    """
    ユーザー明確化の決定スキーマ.
    
    deep_researchの ClarifyWithUser と完全に同じ構造。
    Structured Outputで明確化の要否を判定。
    """
    need_clarification: bool = Field(
        description="追加の明確化が必要かどうか"
    )
    question: str = Field(
        description="ユーザーに尋ねる具体的な質問"
    )
    verification: str = Field(
        description="情報収集完了後の確認メッセージ"
    )

class SimulationBrief(BaseModel):
    """
    シミュレーション要件の構造化スキーマ.
    
    deep_researchの ResearchQuestion に対応。
    会話履歴から構造化されたブリーフを生成。
    """
    # 構造
    pdb_id: Optional[str] = Field(default=None, description="PDB ID")
    fasta_sequence: Optional[str] = Field(default=None, description="FASTA配列")
    ligand_smiles: str = Field(description="リガンドのSMILES")
    
    # パラメータ
    ph: float = Field(default=7.4, description="pH値")
    salt_concentration: float = Field(default=0.15, description="塩濃度 (M)")
    water_model: str = Field(default="TIP3P", description="水モデル")
    box_padding: float = Field(default=12.0, description="Box padding (Å)")
    force_field: str = Field(default="ff19SB", description="力場")
    
    # ワークフロー
    use_boltz2_docking: bool = Field(default=True, description="Boltz-2使用")
    refine_with_smina: bool = Field(default=False, description="Smina精密化")
    output_formats: list[str] = Field(default=["amber"], description="出力形式")
```

```python
# src/mcp_md/state_setup.py (Notebook 2-3で生成)

# ===== セットアップフェーズ状態（Phase 2専用サブグラフ）=====
class SetupState(TypedDict):
    """
    セットアップフェーズのサブグラフ状態.
    
    deep_researchの SupervisorState に対応。
    setup_coordinator と setup_tools ノード間で使用。
    """
    # 入力（親グラフから受け取る）
    simulation_brief: SimulationBrief
    research_brief: str  # deep_researchとの互換性
    
    # 実行状態（サブグラフ内で管理）
    setup_messages: Annotated[Sequence[BaseMessage], add_messages]
    current_step: str
    step_iteration: int  # deep_researchの research_iterations に対応
    
    # 出力（親グラフへ返却）
    outputs: dict
    decision_log: Annotated[list[dict], operator.add]
    raw_notes: Annotated[list[str], operator.add]  # deep_researchと同じ

# ===== セットアップフェーズ出力状態 =====
class SetupOutputState(TypedDict):
    """
    セットアップフェーズのサブグラフ出力.
    
    deep_researchの SupervisorOutputState に対応。
    親グラフに返却するフィールドのみ定義。
    """
    outputs: dict
    decision_log: Annotated[list[dict], operator.add]
    setup_messages: Annotated[Sequence[BaseMessage], add_messages]
    raw_notes: Annotated[list[str], operator.add]

# ===== Structured Tools（Coordinatorパターン用）=====
from langchain_core.tools import tool

@tool
class ExecuteSetupStep(BaseModel):
    """
    固定スケルトンの1ステップを実行するツール.
    
    deep_researchの ConductResearch に対応。
    Coordinator が使用する Structured Tool。
    """
    step_name: str = Field(
        description="実行するステップ名（structure_fetch, ligand_param等）"
    )
    tool_name: str = Field(
        description="使用するMCPツール名"
    )
    parameters: dict = Field(
        description="ツールパラメータ"
    )
    reason: str = Field(
        description="このツール選択の理由（説明可能性）"
    )

@tool
class SetupComplete(BaseModel):
    """
    セットアップ完了を示すツール.
    
    deep_researchの ResearchComplete に対応。
    """
    pass
```

**deep_researchとの対応表**:

| deep_research | mcp-md | 用途 |
|--------------|--------|------|
| `AgentInputState(MessagesState)` | `AgentInputState(MessagesState)` | ユーザー入力 |
| `AgentState(MessagesState)` | `AgentState(MessagesState)` | メイン状態 |
| `SupervisorState(TypedDict)` | `SetupState(TypedDict)` | サブグラフ状態 |
| `SupervisorOutputState` | `SetupOutputState` | サブグラフ出力 |
| `ClarifyWithUser` | `ClarifyWithUser` | 明確化判定 |
| `ResearchQuestion` | `SimulationBrief` | ブリーフ生成 |
| `ConductResearch` | `ExecuteSetupStep` | タスク委譲 |
| `ResearchComplete` | `SetupComplete` | 完了シグナル |

**重要な設計原則**:
1. **入力状態の分離**: `InputState`でユーザー入力のみを受け取る
2. **サブグラフ状態の独立**: サブグラフ専用の`TypedDict`を定義
3. **出力状態の明示**: サブグラフから親への返却フィールドを明示
4. **Structured Output**: Pydanticスキーマで決定を構造化

#### 2. MCP統合の設定

```python
from langchain_mcp_adapters.client import MultiServerMCPClient

# MCPクライアント設定
def create_mcp_client() -> MultiServerMCPClient:
    """MCPクライアントを作成"""
    return MultiServerMCPClient(
        {
            "structure": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.structure_server"]
            },
            "genesis": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.genesis_server"]
            },
            "complex": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.complex_server"]
            },
            "ligand": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.ligand_server"]
            },
            "assembly": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.assembly_server"]
            },
            "export": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.export_server"]
            },
            "qc_min": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.qc_min_server"]
            }
        }
    )

# MCPツールを取得（非同期）
async def load_all_mcp_tools() -> dict[str, Tool]:
    """全MCPツールを読み込み"""
    client = create_mcp_client()
    tools = await client.get_tools()
    # ツール名でアクセス可能なように辞書化
    return {tool.name: tool for tool in tools}
```

**注意**: 
- `MultiServerMCPClient`はデフォルトで**ステートレス**（各ツール呼び出しごとにセッション作成・破棄）
- ステートフルな使用が必要な場合は `async with client.session("server_name")` を使用

#### 3. グラフ構築（3フェーズ統合）

```python
from langchain.chat_models import init_chat_model
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.sqlite import SqliteSaver

# LLMモデル初期化
model = init_chat_model(model="anthropic:claude-sonnet-4-20250514")

# ===== Phase 1: Clarification Subgraph =====
def build_clarification_graph() -> StateGraph:
    """Clarificationフェーズのサブグラフを構築"""
    from core.clarification_nodes import clarify_requirements, generate_simulation_brief
    
    graph = StateGraph(AgentState, input_schema=AgentInputState)
    
    # ノード追加
    graph.add_node("clarify_requirements", clarify_requirements)
    graph.add_node("generate_simulation_brief", generate_simulation_brief)
    
    # エッジ定義（Command APIでルーティング）
    graph.add_edge(START, "clarify_requirements")
    # clarify_requirementsがCommandでルーティング先を決定
    graph.add_edge("generate_simulation_brief", END)
    
    return graph.compile()

# ===== Phase 2: Setup Subgraph =====
async def build_setup_graph() -> StateGraph:
    """Setupフェーズのサブグラフを構築"""
    from core.setup_nodes import setup_coordinator, setup_tools
    
    # MCPツールを読み込み
    mcp_tools = await load_all_mcp_tools()
    
    graph = StateGraph(SetupState, output_schema=SetupOutputState)
    
    # ノード追加（MCPツールを渡す）
    graph.add_node("setup_coordinator", setup_coordinator)
    graph.add_node("setup_tools", setup_tools)
    
    # エッジ定義（Coordinator-Toolsパターン）
    graph.add_edge(START, "setup_coordinator")
    # setup_coordinatorとsetup_toolsがCommandでルーティング
    
    return graph.compile()

# ===== Phase 3: Validation & Export =====
def build_validation_graph() -> StateGraph:
    """Validation & Exportフェーズのサブグラフを構築"""
    from core.validation_nodes import validate_system, export_files, generate_report
    
    graph = StateGraph(AgentState)
    
    # ノード追加
    graph.add_node("validate_system", validate_system)
    graph.add_node("export_files", export_files)
    graph.add_node("generate_report", generate_report)
    
    # エッジ定義（直線的フロー）
    graph.add_edge(START, "validate_system")
    graph.add_edge("validate_system", "export_files")
    graph.add_edge("export_files", "generate_report")
    graph.add_edge("generate_report", END)
    
    return graph.compile()

# ===== メイングラフ（3フェーズ統合）=====
async def create_agent() -> StateGraph:
    """メインエージェントグラフを構築"""
    
    # サブグラフを構築
    clarification_graph = build_clarification_graph()
    setup_graph = await build_setup_graph()
    validation_graph = build_validation_graph()
    
    # メイングラフ
    main_graph = StateGraph(AgentState, input_schema=AgentInputState)
    
    # サブグラフをノードとして追加
    main_graph.add_node("clarification_phase", clarification_graph)
    main_graph.add_node("setup_phase", setup_graph)
    main_graph.add_node("validation_phase", validation_graph)
    
    # フェーズ間のエッジ定義
    main_graph.add_edge(START, "clarification_phase")
    main_graph.add_edge("clarification_phase", "setup_phase")
    main_graph.add_edge("setup_phase", "validation_phase")
    main_graph.add_edge("validation_phase", END)
    
    # チェックポイント機能（永続化）
    memory = SqliteSaver.from_conn_string("checkpoints/workflow.db")
    
    # コンパイル
    return main_graph.compile(
        checkpointer=memory,
        interrupt_before=["setup_phase"]  # セットアップ前に人間確認（オプション）
    )

# エージェント作成
agent = await create_agent()
```

#### 4. エージェント実行例

```python
from langchain_core.messages import HumanMessage

# エージェント作成
agent = await create_agent()

# スレッドID（会話セッション管理）
config = {"configurable": {"thread_id": "session_20250103_001"}}

# ===== 実行例1: シンプルな実行 =====
result = await agent.ainvoke(
    {
        "messages": [HumanMessage(content="PDB 1ABCにAspirinをドッキングしてAmber形式で出力して")]
    },
    config=config
)

print(result["final_report"])

# ===== 実行例2: 対話的な明確化 =====
# 最初の入力（情報不足）
result1 = await agent.ainvoke(
    {
        "messages": [HumanMessage(content="タンパク質にリガンドをドッキングしたい")]
    },
    config=config
)
# AIからの質問: "どのタンパク質ですか？PDB IDまたはFASTA配列を教えてください"

# ユーザーの追加情報
result2 = await agent.ainvoke(
    {
        "messages": [HumanMessage(content="PDB ID は 7BV2 で、リガンドはCC(=O)Oc1ccccc1C(=O)O")]
    },
    config=config
)
# AIが自動でセットアップを開始

# ===== 実行例3: 中断と再開 =====
# Setup前で中断（interrupt_before=["setup_phase"]）
result = await agent.ainvoke(
    {
        "messages": [HumanMessage(content="PDB 1ABC, Aspirin")]
    },
    config=config
)

# 現在の状態を確認
current_state = agent.get_state(config)
print(f"Simulation Brief: {current_state.values['simulation_brief']}")

# ユーザー承認後、再開
user_approval = input("Proceed with setup? (y/n): ")
if user_approval == "y":
    result = await agent.ainvoke(None, config=config)

# ===== 実行例4: 状態履歴の確認 =====
# 過去の実行履歴を取得（time-travel）
for state in agent.get_state_history(config):
    print(f"Checkpoint ID: {state.config['configurable']['checkpoint_id']}")
    print(f"Current Step: {state.values.get('current_step', 'N/A')}")
    print(f"Decision Log: {len(state.values.get('decision_log', []))} decisions")
    print("---")

# 特定のチェックポイントから再開
checkpoint_config = {
    "configurable": {
        "thread_id": "session_20250103_001",
        "checkpoint_id": "specific_checkpoint_id"
    }
}
result = await agent.ainvoke(None, config=checkpoint_config)
```

#### 5. ストリーミング実行（リアルタイム進捗表示）

```python
# ストリーミングモードで実行
async for event in agent.astream(
    {
        "messages": [HumanMessage(content="PDB 1ABC, Aspirin")]
    },
    config=config,
    stream_mode="updates"  # "values", "updates", "messages"から選択
):
    # 各ノードの更新をリアルタイム表示
    for node_name, node_update in event.items():
        if node_name == "setup_coordinator":
            print(f"🤖 Coordinator: {node_update.get('setup_messages', [])[-1].content}")
        elif node_name == "setup_tools":
            print(f"🔧 Tool executed: {node_update.get('decision_log', [])[-1]}")
```

#### 6. 安全性・認証
- **MCP統合**: MCPサーバーは標準プロトコルで動作（認証は別管理）
- **実装方針**:
  - 環境変数で機密値（API key等）を管理
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

## 6. ロードマップ（段階的実装）

### 実装ステップ（deep_research_from_scratchパターン完全適用）

#### Step 1: Notebookベース開発（4週間）🎯

**目標**: 各フェーズをJupyter Notebookで段階的に実装・検証

**重要な開発原則**（deep_researchから学ぶ）:
```
🚨 CRITICAL: Notebooks are the source of truth!

notebooks/     ← MODIFY THESE (開発の中心) ✏️
src/mcp_md/    ← GENERATED CODE (直接編集禁止) 🚫
```

**Notebook構成**（5段階・deep_researchパターン）:

##### 1. `1_clarification.ipynb` - Phase 1: Scoping
**deep_researchの`1_scoping.ipynb`に対応**

実装内容：
- **Structured Output定義**:
  ```python
  class ClarifyWithUser(BaseModel):
      need_clarification: bool
      question: str
      verification: str
  
  class SimulationBrief(BaseModel):
      pdb_id: Optional[str]
      fasta_sequence: Optional[str]
      ligand_smiles: str
      ph: float = 7.4
      # ... その他パラメータ
  ```

- **2ノードワークフロー**:
  ```python
  def clarify_requirements(state: AgentState) -> Command[...]:
      # Structured Outputで明確化判定
      model_with_structured_output = model.with_structured_output(ClarifyWithUser)
      response = model_with_structured_output.invoke(...)
      
      if response.need_clarification:
          return Command(goto=END, update={"messages": [...]})
      else:
          return Command(goto="generate_simulation_brief", update={...})
  
  def generate_simulation_brief(state: AgentState):
      # Structured Outputでブリーフ生成
      model_with_structured_output = model.with_structured_output(SimulationBrief)
      brief = model_with_structured_output.invoke(...)
      return {"simulation_brief": brief, "setup_messages": [...]}
  ```

- **`%%writefile`で以下を生成**:
  - `src/mcp_md/state_scope.py` - 状態定義
  - `src/mcp_md/clarification_agent.py` - Clarificationエージェント
  - `src/mcp_md/prompts.py` - プロンプトテンプレート（段階的に追加）

**deep_researchとの対応**:
| deep_research | mcp-md | 目的 |
|--------------|--------|------|
| `clarify_with_user` | `clarify_requirements` | 明確化判定 |
| `write_research_brief` | `generate_simulation_brief` | ブリーフ生成 |
| `ResearchQuestion` | `SimulationBrief` | 構造化ブリーフ |

---

##### 2. `2_setup_agent.ipynb` - Phase 2: Setup (基本実装)
**deep_researchの`2_research_agent.ipynb`に対応**

実装内容：
- **MCPツール統合**（7サーバー）:
  ```python
  from langchain_mcp_adapters.client import MultiServerMCPClient
  
  async def load_mcp_tools():
      client = MultiServerMCPClient({
          "structure": {"transport": "stdio", "command": "python", "args": ["-m", "servers.structure_server"]},
          "genesis": {...},
          # ... 他5サーバー
      })
      tools = await client.get_tools()
      return {tool.name: tool for tool in tools}
  ```

- **固定スケルトン実装**（シンプルな直線的フロー）:
  ```python
  # 6ステップの固定スケルトン
  SETUP_STEPS = [
      "structure_fetch",
      "structure_repair", 
      "ligand_param",
      "complex_generation",
      "assembly",
      "qc_check"
  ]
  
  # シンプルなReActエージェント
  def setup_agent(state: SetupState):
      # 現在のステップに応じたツール選択
      step = state["current_step"]
      available_tools = get_tools_for_step(step, mcp_tools)
      
      model_with_tools = model.bind_tools(available_tools)
      response = model_with_tools.invoke(state["setup_messages"])
      
      return {"setup_messages": [response]}
  
  def setup_tools(state: SetupState):
      # ツール実行（同期的・シンプル）
      results = []
      for tool_call in state["setup_messages"][-1].tool_calls:
          result = mcp_tools[tool_call["name"]].invoke(tool_call["args"])
          results.append(ToolMessage(content=result, ...))
      
      # 次のステップへ
      next_step = get_next_step(state["current_step"])
      return {"setup_messages": results, "current_step": next_step}
  ```

- **`%%writefile`で以下を生成**:
  - `src/mcp_md/state_setup.py` - Setup状態定義
  - `src/mcp_md/setup_agent.py` - 基本Setupエージェント
  - `src/mcp_md/mcp_integration.py` - MCP統合
  - `src/mcp_md/utils.py` - ユーティリティ関数

**deep_researchとの対応**:
| deep_research | mcp-md | 目的 |
|--------------|--------|------|
| `researcher` | `setup_agent` | エージェントノード |
| `researcher_tools` | `setup_tools` | ツール実行ノード |
| Tavily search | MCP servers (7つ) | 外部ツール |

---

##### 3. `3_setup_coordinator.ipynb` - Phase 2: Setup (Coordinator-Tools高度化)
**deep_researchの`4_research_supervisor.ipynb`に対応**

実装内容：
- **Coordinator-Toolsパターン**（Supervisorパターン適用）:
  ```python
  # Structured Tools定義
  @tool
  class ExecuteSetupStep(BaseModel):
      """固定スケルトンの1ステップを実行"""
      step_name: str = Field(description="実行するステップ名")
      tool_name: str = Field(description="使用するMCPツール名")
      parameters: dict = Field(description="ツールパラメータ")
      reason: str = Field(description="このツール選択の理由")
  
  @tool
  class SetupComplete(BaseModel):
      """セットアップ完了を示す"""
      pass
  
  # Coordinator ノード（決定）
  async def setup_coordinator(state: SetupState) -> Command[Literal["setup_tools"]]:
      """固定スケルトンに従い、現在のステップに最適なツールを選択"""
      setup_tools = [ExecuteSetupStep, SetupComplete, think_tool]
      model_with_tools = model.bind_tools(setup_tools)
      
      system_prompt = setup_coordinator_prompt.format(
          current_step=state["current_step"],
          simulation_brief=state["simulation_brief"],
          available_tools=get_tools_for_step(state["current_step"])
      )
      
      response = await model_with_tools.ainvoke(
          [SystemMessage(content=system_prompt)] + state["setup_messages"]
      )
      
      return Command(
          goto="setup_tools",
          update={"setup_messages": [response]}
      )
  
  # Tools ノード（実行）
  async def setup_tools(state: SetupState) -> Command[Literal["setup_coordinator", "__end__"]]:
      """Coordinatorの決定を実行し、結果を記録"""
      most_recent_message = state["setup_messages"][-1]
      
      if not most_recent_message.tool_calls:
          return Command(goto="setup_coordinator")
      
      tool_results = []
      decision_logs = []
      
      for tool_call in most_recent_message.tool_calls:
          if tool_call["name"] == "think_tool":
              # think_tool実行
              result = think_tool.invoke(tool_call["args"])
              tool_results.append(ToolMessage(content=result, ...))
              
          elif tool_call["name"] == "ExecuteSetupStep":
              # MCPツール実行
              mcp_tool_name = tool_call["args"]["tool_name"]
              result = await mcp_tools[mcp_tool_name].ainvoke(tool_call["args"]["parameters"])
              tool_results.append(ToolMessage(content=result, ...))
              
              # 決定ログ記録
              decision_logs.append({
                  "step": tool_call["args"]["step_name"],
                  "tool": mcp_tool_name,
                  "parameters": tool_call["args"]["parameters"],
                  "reason": tool_call["args"]["reason"],
                  "timestamp": datetime.now().isoformat()
              })
              
          elif tool_call["name"] == "SetupComplete":
              # セットアップ完了
              return Command(
                  goto=END,
                  update={
                      "setup_messages": tool_results,
                      "decision_log": decision_logs
                  }
              )
      
      # 次のステップへ
      next_step = advance_step(state["current_step"], decision_logs)
      return Command(
          goto="setup_coordinator",
          update={
              "setup_messages": tool_results,
              "decision_log": decision_logs,
              "current_step": next_step
          }
      )
  ```

- **`%%writefile`で以下を更新**:
  - `src/mcp_md/setup_coordinator.py` - Coordinatorエージェント
  - `src/mcp_md/prompts.py` - Coordinatorプロンプト追加
  - `src/mcp_md/decision_logger.py` - 決定ログ管理

**deep_researchとの対応**:
| deep_research | mcp-md | 目的 |
|--------------|--------|------|
| `supervisor` | `setup_coordinator` | タスク調整・ツール選択 |
| `supervisor_tools` | `setup_tools` | ツール実行・結果集約 |
| `ConductResearch` | `ExecuteSetupStep` | タスク委譲ツール |
| `ResearchComplete` | `SetupComplete` | 完了シグナル |
| 並列research agents | 順次setup steps | 実行モデル（違いに注意） |

**重要な違い**:
- deep_research: 並列研究（`asyncio.gather()`で複数のresearch agentを同時実行）
- mcp-md: **順次実行**（固定スケルトンの各ステップを順番に実行）

---

##### 4. `4_validation.ipynb` - Phase 3: Validation & Export
**deep_researchの`5_full_agent.ipynb`のWrite部分に対応**

実装内容：
- **3ノードの直線的フロー**:
  ```python
  def validate_system(state: AgentState):
      """QC検証実行"""
      qc_results = await run_full_qc(
          prmtop=state["outputs"]["prmtop"],
          inpcrd=state["outputs"]["inpcrd"]
      )
      return {"qc_results": qc_results, ...}
  
  def export_files(state: AgentState):
      """形式変換・エクスポート"""
      exports = {}
      for format in state["simulation_brief"]["output_formats"]:
          exports[format] = await convert_to_format(state["outputs"], format)
      return {"exports": exports, ...}
  
  def generate_report(state: AgentState):
      """最終レポート生成"""
      report = {
          "simulation_brief": state["simulation_brief"].model_dump(),
          "decision_log": state["decision_log"],
          "qc_results": state["qc_results"],
          ...
      }
      markdown_report = await generate_markdown_report(report)
      return {"final_report": markdown_report, ...}
  ```

- **`%%writefile`で以下を生成**:
  - `src/mcp_md/validation_agent.py` - Validationエージェント
  - `src/mcp_md/report_generation.py` - レポート生成
  - `src/mcp_md/prompts.py` - レポート生成プロンプト追加

---

##### 5. `5_full_agent.ipynb` - 全統合（End-to-End）
**deep_researchの`5_full_agent.ipynb`に完全対応**

実装内容：
- **3フェーズサブグラフ統合**:
  ```python
  # 各フェーズをサブグラフとして構築
  clarification_graph = build_clarification_graph()
  setup_graph = await build_setup_graph()
  validation_graph = build_validation_graph()
  
  # メイングラフ
  main_graph = StateGraph(AgentState, input_schema=AgentInputState)
  
  # サブグラフをノードとして追加
  main_graph.add_node("clarification_phase", clarification_graph)
  main_graph.add_node("setup_phase", setup_graph)
  main_graph.add_node("validation_phase", validation_graph)
  
  # フェーズ間のエッジ
  main_graph.add_edge(START, "clarification_phase")
  main_graph.add_edge("clarification_phase", "setup_phase")
  main_graph.add_edge("setup_phase", "validation_phase")
  main_graph.add_edge("validation_phase", END)
  
  # チェックポイント
  from langgraph.checkpoint.sqlite import SqliteSaver
  memory = SqliteSaver.from_conn_string("checkpoints/workflow.db")
  
  # コンパイル
  agent = main_graph.compile(
      checkpointer=memory,
      interrupt_before=["setup_phase"]  # オプション: 人間確認
  )
  ```

- **対話的実行**:
  ```python
  # スレッドベースのセッション管理
  config = {"configurable": {"thread_id": "session_20250103_001"}}
  
  # 1. 初回実行（情報不足で明確化）
  result1 = await agent.ainvoke({
      "messages": [HumanMessage(content="タンパク質にリガンドをドッキングしたい")]
  }, config)
  # → AI: "どのタンパク質ですか？PDB IDまたはFASTA配列を教えてください"
  
  # 2. 追加情報で再実行
  result2 = await agent.ainvoke({
      "messages": [HumanMessage(content="PDB 7BV2で、リガンドはアスピリン")]
  }, config)
  # → AI: 自動でセットアップ開始
  
  # 3. ストリーミングで進捗確認
  async for event in agent.astream({...}, config, stream_mode="updates"):
      for node_name, node_update in event.items():
          print(f"[{node_name}] {node_update}")
  ```

- **`%%writefile`で以下を生成**:
  - `src/mcp_md/full_agent.py` - 統合エージェント
  - `src/mcp_md/__init__.py` - パッケージエクスポート

---

**開発ワークフロー**（deep_research完全準拠）:
```
1. Notebookで開発 (notebooks/*.ipynb)
   ├─ インタラクティブに実装・テスト
   ├─ %%writefile でsrc/mcp_md/にコード生成
   └─ 実行結果をその場で確認

2. 生成されたコードを利用 (src/mcp_md/)
   ├─ パッケージとしてインポート可能
   ├─ 本番実行用
   └─ 直接編集禁止🚫

3. 修正が必要な場合
   ├─ Notebookに戻る ✅
   ├─ %%writefileセルを修正
   └─ セル再実行でsrc/を再生成
```

**成果物**:
- ✅ 動作するNotebook群（5個）
- ✅ 自動生成されたソースコード（`src/mcp_md/`）
- ✅ 実装ドキュメント
- ✅ 実行履歴とデモ結果

#### Step 2: ソースコードパッケージ化（2週間）

**目標**: Notebookから生成されたコードをパッケージ化

**実装内容**:
```
src/mcp_md/
├── __init__.py
├── states.py                    # 状態定義
├── schemas.py                   # Structured Output スキーマ
├── prompts.py                   # プロンプトテンプレート
├── utils.py                     # ユーティリティ
├── clarification_agent.py       # Phase 1
├── setup_agent.py               # Phase 2
├── validation_agent.py          # Phase 3
└── full_agent.py                # 統合エージェント
```

**機能**:
- Amber特化パイプライン（PDB/FASTA → prmtop/inpcrd）
- 基本QC（最小化、電荷整合）
- MCPサーバー統合（7サーバー）

#### Step 3: QC強化（2-3週間）

**目標**: 学術発表レベルの品質保証

**追加機能**:
- 自作MolProbityによる物理化学的一貫性チェック
- PoseBustersチェック統合
- QCレポート生成（JSON + Markdown）
- QC失敗時の自動リトライ

#### Step 4: 出力拡張（3-4週間）

**目標**: GROMACS/OpenMM対応、膜系

**追加機能**:
- GROMACS出力（ParmEd）
- OpenMM XML出力
- 膜系構築（Packmol-Memgen）
- 混合溶媒対応（Packmol）

#### Step 5: HPC/永続運用（4-6週間）

**目標**: 大規模スクリーニング、長期運用

**追加機能**:
- LangGraph Studio統合（可視化）
- 結果キャッシュ（同一入力の再利用）
- 結果索引（検索可能なDB）
- HPC連携（Slurmジョブ投入）

#### Step 6: 将来拡張（低優先度）

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

## 7. 現在の実装状況

### 実装済み（7 FastMCP Servers）✅

| Component | Status | 主要機能 |
|-----------|--------|---------|
| Structure Server | ✅ | PDB取得、PDBFixer、PDB2PQR、構造検証 |
| Genesis Server | ✅ | Boltz-2タンパク質生成（FASTA→PDB、マルチマー） |
| Complex Server | ✅ | Boltz-2複合体予測、Sminaドッキング、ポーズ精密化 |
| Ligand Server | ✅ | RDKit 3D生成、AmberTools GAFF2/AM1-BCC |
| Assembly Server | ✅ | tleap系構築、Packmol-Memgen膜系 |
| Export Server | ✅ | Amber/GROMACS/OpenMM形式変換、パッケージング |
| QC/Min Server | ✅ | OpenMM最小化、衝突検出、結合長・キラリティチェック |

### 未実装（LangGraph Agent）❌

| Component | Status | 実装予定 |
|-----------|--------|---------|
| Clarification Phase | ❌ | Step 1 (Notebook 1) |
| Setup Phase | ❌ | Step 1 (Notebook 2-3) |
| Validation Phase | ❌ | Step 1 (Notebook 4) |
| Full Agent Integration | ❌ | Step 1 (Notebook 5) |
| Structured Output Schemas | ❌ | Step 1 (全Notebook) |
| Command API Integration | ❌ | Step 1 (全Notebook) |
| Checkpoint/Persistence | ❌ | Step 1 (Notebook 5) |
| Prompts Library | ❌ | Step 1 (全Notebook) |

### 新しいディレクトリ構造（deep_researchパターン完全適用）

#### 推奨構造（Notebookベース開発）

**重要**: deep_research_from_scratchと全く同じワークフローを採用

```
mcp-md/
│
├── notebooks/                    # 🎯 開発の中心（MODIFY THESE）✏️
│   │                             # deep_researchの5段階構成を完全に踏襲
│   ├── 1_clarification.ipynb     # Phase 1: Scoping
│   │                             # └→ clarify_requirements + generate_simulation_brief
│   ├── 2_setup_agent.ipynb       # Phase 2: Setup (基本・ReAct)
│   │                             # └→ setup_agent + setup_tools (シンプル)
│   ├── 3_setup_coordinator.ipynb # Phase 2: Setup (Coordinator-Tools)
│   │                             # └→ setup_coordinator + setup_tools (高度)
│   ├── 4_validation.ipynb        # Phase 3: Validation & Export
│   │                             # └→ validate + export + report
│   ├── 5_full_agent.ipynb        # 全統合 (End-to-End)
│   │                             # └→ 3フェーズサブグラフ統合
│   └── utils.py                  # Notebook用ユーティリティ（rich表示等）
│
├── src/mcp_md/                   # 🚫 生成されたソースコード（DO NOT MODIFY）
│   │                             # %%writefileで自動生成される
│   ├── __init__.py               # パッケージエクスポート
│   │
│   # 状態定義（deep_researchパターン）
│   ├── state_scope.py            # AgentInputState, AgentState, ClarifyWithUser, SimulationBrief
│   ├── state_setup.py            # SetupState, SetupOutputState
│   │
│   # エージェント実装（段階的に生成）
│   ├── clarification_agent.py    # Phase 1: clarify_requirements + generate_simulation_brief
│   ├── setup_agent.py            # Phase 2基本: setup_agent + setup_tools (ReAct)
│   ├── setup_coordinator.py      # Phase 2高度: setup_coordinator + setup_tools (Supervisor)
│   ├── validation_agent.py       # Phase 3: validate + export + report
│   ├── full_agent.py             # 統合: 3フェーズサブグラフ
│   │
│   # 共通モジュール
│   ├── prompts.py                # プロンプトテンプレート（段階的に追加）
│   ├── mcp_integration.py        # MCP統合（MultiServerMCPClient）
│   ├── decision_logger.py        # 決定ログ管理
│   ├── report_generation.py      # レポート生成
│   └── utils.py                  # ユーティリティ関数
│
├── servers/                      # ✅ FastMCPサーバー（既存・維持）
│   ├── __init__.py
│   ├── structure_server.py       # PDB取得、修復、プロトネーション
│   ├── genesis_server.py         # Boltz-2タンパク質生成
│   ├── complex_server.py         # Boltz-2複合体、Sminaドッキング
│   ├── ligand_server.py          # GAFF2/AM1-BCC パラメータ化
│   ├── assembly_server.py        # tleap系構築、膜系
│   ├── export_server.py          # Amber/GROMACS/OpenMM形式変換
│   └── qc_min_server.py          # 最小化、QC検証
│
├── common/                       # ✅ 共通ライブラリ（既存・維持）
│   ├── __init__.py
│   ├── base.py                   # BaseToolWrapper
│   └── utils.py                  # ロガー、ディレクトリ管理
│
├── checkpoints/                  # LangGraphチェックポイント
│   └── workflow.db               # SQLiteベースの永続化
│
├── runs/                         # 実行結果（タイムスタンプ単位）
│   └── <timestamp>/
│       ├── simulation_brief.json # Phase 1出力
│       ├── decision_log.json     # Phase 2決定ログ
│       ├── outputs/              # PDB, prmtop, inpcrd等
│       ├── qc_report.json        # Phase 3検証結果
│       └── metadata.json         # 再現性情報
│
├── pyproject.toml                # 依存関係（uv管理）
├── README.md                     # プロジェクト説明
├── ARCHITECTURE.md               # このファイル
├── AGENTS.md                     # Cursor AI Agent設定
└── .cursor/                      # Cursorプロジェクト設定
    └── rules/                    # プロジェクトルール
        ├── project-rules.md      # プロジェクト全体のルール
        └── notebook-development.md # Notebook開発ルール
```

#### 開発ワークフロー（deep_research完全準拠）

**🚨 CRITICAL: ファイル編集のルール**

```
✅ DO:     notebooks/*.ipynb を編集
✅ DO:     %%writefile で src/mcp_md/ を生成
✅ DO:     Notebookでテスト・実行

🚫 DON'T: src/mcp_md/ を直接編集
🚫 DON'T: 手動で src/ ファイルを作成
```

**開発サイクル**:

1. **Notebookで実装** (`notebooks/*.ipynb`)
   ```python
   # Notebookセル
   %%writefile ../src/mcp_md/clarification_agent.py
   
   """Clarification Agent Implementation"""
   
   from langchain.chat_models import init_chat_model
   # ... 実装コード ...
   ```

2. **その場でテスト** (同じNotebook内)
   ```python
   # 次のセルで即座にテスト
   from mcp_md.clarification_agent import clarify_requirements
   
   result = await clarify_requirements(test_state)
   print(result)
   ```

3. **修正が必要な場合**
   - ❌ `src/mcp_md/clarification_agent.py`を直接編集 → ダメ！
   - ✅ Notebookの`%%writefile`セルに戻る → 修正 → 再実行

**メリット**:
- インタラクティブな開発（即座にテスト）
- バージョン管理がシンプル（Notebookが真実の源泉）
- ドキュメント統合（説明とコードが一体）
- 再現性（Notebookを実行すればsrc/が生成される）

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

#### 3. LangGraph × MCP統合パターン

```python
# core/mcp_integration.py
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_core.tools import Tool

def create_mcp_client() -> MultiServerMCPClient:
    """MCPクライアントを作成"""
    return MultiServerMCPClient(
        {
            "structure": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.structure_server"]
            },
            "genesis": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.genesis_server"]
            },
            "complex": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.complex_server"]
            },
            "ligand": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.ligand_server"]
            },
            "assembly": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.assembly_server"]
            },
            "export": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.export_server"]
            },
            "qc_min": {
                "transport": "stdio",
                "command": "python",
                "args": ["-m", "servers.qc_min_server"]
            }
        }
    )

async def load_all_mcp_tools() -> dict[str, Tool]:
    """全MCPツールを読み込み"""
    client = create_mcp_client()
    tools = await client.get_tools()
    # ツール名でアクセス可能なように辞書化
    return {tool.name: tool for tool in tools}

# ステートフルな使用が必要な場合
async def load_mcp_tools_stateful(server_name: str):
    """特定のサーバーからステートフルにツールを読み込み"""
    from langchain_mcp_adapters.tools import load_mcp_tools
    
    client = create_mcp_client()
    async with client.session(server_name) as session:
        tools = await load_mcp_tools(session)
        return tools

# core/workflow_graph.py
from langgraph.graph import StateGraph, END
from langgraph.checkpoint.sqlite import SqliteSaver
from .workflow_state import WorkflowState
from .workflow_nodes import (
    planner_node,
    create_structure_fetch_node,
    create_repair_node,
    # ...
)
from .mcp_integration import load_all_mcp_tools

async def create_workflow_graph():
    """ワークフローグラフを構築"""
    # MCPツール読み込み
    mcp_tools = await load_all_mcp_tools()
    
    # グラフ構築
    graph = StateGraph(WorkflowState)
    
    # ノード追加（ツールを渡す）
    graph.add_node("planner", planner_node)
    graph.add_node("fetch", create_structure_fetch_node(mcp_tools))
    graph.add_node("repair", create_repair_node(mcp_tools))
    # ...
    
    # エッジ定義
    graph.set_entry_point("planner")
    graph.add_edge("planner", "fetch")
    # ...
    graph.add_edge("qc", END)
    
    # チェックポイント設定
    memory = SqliteSaver.from_conn_string("checkpoints/workflow.db")
    
    return graph.compile(checkpointer=memory)
```

**MCPトランスポートタイプ**:
- **stdio**: ローカルサブプロセス通信（本プロジェクトで使用）
- **streamable_http**: HTTPベースのリモートサーバー
- **SSE (Server-Sent Events)**: リアルタイムストリーミング通信用

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

### 環境セットアップ（推奨: conda + uv）

#### 1. conda環境作成と科学計算ツールのインストール

```bash
# conda環境作成
conda create -n mcp-md python=3.11
conda activate mcp-md

# 科学計算ツール（conda-forge推奨）
conda install -c conda-forge ambertools packmol smina pdbfixer

# uvのインストール（conda環境内）
pip install uv
```

#### 2. conda環境内でuvを使ってPythonパッケージをインストール

```bash
# conda環境がアクティブな状態で実行
conda activate mcp-md

# 基本依存関係のインストール（conda環境に直接インストール）
uv pip install -e .

# または、pyproject.tomlから直接インストール
uv pip install --project pyproject.toml

# 特定のLLMプロバイダーも含める場合
uv pip install -e ".[openai]"      # OpenAI/LM Studio
uv pip install -e ".[anthropic]"   # Claude
uv pip install -e ".[google]"      # Gemini

# 全てのオプション依存関係
uv pip install -e ".[openai,anthropic,google,dev]"
```

#### 3. 実行方法

```bash
# conda環境がアクティブな状態で
conda activate mcp-md

# uv runを使って実行（高速起動）
uv run python main.py

# またはMCPサーバーの起動
uv run python -m servers.structure_server

# LangGraphワークフローの実行
uv run python -m core.workflow_graph

# 通常のpythonコマンドも使用可能
python main.py
```

#### 4. pyproject.toml設定例

```toml
[project]
name = "mcp-md"
version = "0.1.0"
description = "Amber-focused MD setup with LangGraph + MCP"
requires-python = ">=3.11"
dependencies = [
    "boltz>=2.0.0",
    "pdb2pqr>=3.1.0",
    "propka>=3.5.0",
    "rdkit>=2023.9.1",
    "openmm>=8.3.1",
    "parmed>=4.3.0",
    "fastmcp>=0.1.0",
    "langchain-core>=1.0.0",
    "langgraph>=0.2.0",
    "langchain-mcp-adapters>=0.1.0",  # MCP統合
]

[project.optional-dependencies]
openai = ["langchain-openai>=0.2.0"]
anthropic = ["langchain-anthropic>=0.3.0"]
google = ["langchain-google-genai>=0.1.0"]
dev = ["pytest>=7.0", "black>=24.0", "ruff>=0.1.0"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
```

### 主要パッケージ一覧

#### 科学計算ツール（conda経由）
- **AmberTools**: 完全OSSのAmberツール群（tleap, antechamber, parmchk2等）
- **Packmol**: 溶媒・膜系の構築
- **Smina**: ドッキングツール（AutoDock Vina fork）
- **PDBFixer**: PDB構造修復

#### Pythonパッケージ（uv経由）
- **Boltz-2**: 構造予測・複合体生成
- **PDB2PQR + PROPKA**: プロトネーション
- **RDKit**: ケモインフォマティクス
- **OpenMM + ParmEd**: MD計算とトポロジー変換
- **FastMCP**: MCPサーバー実装
- **LangChain Core + LangGraph**: ワークフロー構築

### 注意事項

1. **conda + uv併用の方針**: 
   - **conda環境**: 科学計算ツール（C/C++バイナリ）+ Python本体
   - **uv pip**: conda環境内でPythonパッケージをインストール（高速）
   - **uv run**: conda環境内でスクリプト実行（キャッシュ活用で高速起動）
   
2. **uv独自の仮想環境は使わない**: 
   - `uv sync` は実行しない（`.venv`を作成してしまう）
   - `uv pip install` を使ってconda環境に直接インストール
   - `uv run` はconda環境のPythonを使用
   
3. **依存関係のロック**: 
   - conda環境では `conda env export > environment.yml` でロック
   - Pythonパッケージは `uv pip compile pyproject.toml -o requirements.txt` でロック可能
   - または `pip freeze > requirements.txt`

4. **MCP統合**: 
   - `langchain-mcp-adapters`パッケージを使用（公式サポート）
   - `MultiServerMCPClient`で複数のMCPサーバーを統合
   - デフォルトはステートレス、必要に応じて`client.session()`でステートフル化

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

#### 主要フレームワーク
- **LangChain**: https://github.com/langchain-ai/langchain
  - **LangChain 1.0 ドキュメント**: https://python.langchain.com/docs/
  - **LangChain 1.0 移行ガイド**: https://python.langchain.com/docs/versions/v0_3/migrating_chains/
  - **MCP統合ドキュメント**: https://docs.langchain.com/oss/python/langchain/mcp
- **LangGraph**: https://github.com/langchain-ai/langgraph
  - **LangGraph ドキュメント**: https://langchain-ai.github.io/langgraph/
  - **チェックポイント機能**: https://langchain-ai.github.io/langgraph/concepts/persistence/
- **langchain-mcp-adapters**: LangChainとMCPの公式統合パッケージ
- **FastMCP**: https://github.com/jlowin/fastmcp
- **MCP Protocol**: https://modelcontextprotocol.io/
- **uv**: https://github.com/astral-sh/uv
  - **uvドキュメント**: https://docs.astral.sh/uv/

#### 科学計算ツール
- **Boltz-2**: https://github.com/jwohlwend/boltz
- **AmberTools**: https://ambermd.org/AmberTools.php
- **OpenMM**: https://openmm.org/
- **PoseBusters**: https://github.com/maabuu/posebusters

---

## 11. まとめ

### 現在地（2025年11月3日時点）

#### 実装済み ✅
- ✅ **7つのFastMCPサーバー**（基本機能完成）
  - Structure, Genesis, Complex, Ligand, Assembly, Export, QC/Min
- ✅ **共通ライブラリ**（`common/`）完成
  - BaseToolWrapper, ロガー、ユーティリティ
- ✅ **deep_research_from_scratch**の参照実装を入手・分析完了

#### 未実装 ❌
- ❌ **Notebookベース開発環境**（最優先）
  - 5つのNotebook（1_clarification → 5_full_agent）
  - notebooks/utils.py（rich表示等）
- ❌ **LangGraphエージェント統合**
  - 状態定義（state_scope.py, state_setup.py）
  - 各フェーズの実装（3フェーズ）
- ❌ **プロンプトライブラリ**
  - Clarification, Setup Coordinator, Validation用

### 次のステップ（deep_researchパターン完全適用）

#### フェーズ1: Notebook環境構築（1週間）
```bash
# 1. notebooks/ディレクトリ作成
mkdir -p notebooks

# 2. deep_researchのutils.pyを参考に作成
# notebooks/utils.py (rich formatting)

# 3. 依存関係追加
uv add jupyter rich langchain-core langgraph langchain-mcp-adapters

# 4. Jupyter起動
uv run jupyter notebook
```

**成果物**:
- notebooks/utils.py（rich表示）
- 実行可能なJupyter環境

---

#### フェーズ2: Notebook 1実装（3-4日）
**1_clarification.ipynb**

実装内容：
1. Structured Output定義（ClarifyWithUser, SimulationBrief）
2. clarify_requirements ノード（Command APIでルーティング）
3. generate_simulation_brief ノード
4. サブグラフ構築とテスト

**`%%writefile`で生成**:
- src/mcp_md/state_scope.py
- src/mcp_md/clarification_agent.py
- src/mcp_md/prompts.py（初版）

---

#### フェーズ3: Notebook 2実装（3-4日）
**2_setup_agent.ipynb**

実装内容：
1. MCP統合（MultiServerMCPClient）
2. 固定スケルトン実装（6ステップ）
3. シンプルなReActエージェント（setup_agent + setup_tools）
4. サブグラフ構築とテスト

**`%%writefile`で生成**:
- src/mcp_md/state_setup.py
- src/mcp_md/setup_agent.py
- src/mcp_md/mcp_integration.py
- src/mcp_md/utils.py

---

#### フェーズ4: Notebook 3実装（4-5日）
**3_setup_coordinator.ipynb**

実装内容：
1. Coordinator-Toolsパターン（Supervisor適用）
2. Structured Tools（ExecuteSetupStep, SetupComplete）
3. 決定ログ記録
4. think_tool統合

**`%%writefile`で更新**:
- src/mcp_md/setup_coordinator.py
- src/mcp_md/prompts.py（Coordinatorプロンプト追加）
- src/mcp_md/decision_logger.py

---

#### フェーズ5: Notebook 4実装（2-3日）
**4_validation.ipynb**

実装内容：
1. QC検証ノード
2. 形式変換・エクスポートノード
3. レポート生成ノード
4. サブグラフ構築とテスト

**`%%writefile`で生成**:
- src/mcp_md/validation_agent.py
- src/mcp_md/report_generation.py
- src/mcp_md/prompts.py（レポート生成プロンプト追加）

---

#### フェーズ6: Notebook 5実装（3-4日）
**5_full_agent.ipynb**

実装内容：
1. 3フェーズサブグラフ統合
2. チェックポイント機能（SqliteSaver）
3. スレッドベースのセッション管理
4. エンドツーエンドテスト
5. ストリーミング実行

**`%%writefile`で生成**:
- src/mcp_md/full_agent.py
- src/mcp_md/__init__.py

---

#### マイルストーン

| フェーズ | 期間 | 成果物 |
|--------|------|--------|
| 1. 環境構築 | 1週間 | Jupyter + utils.py |
| 2. Notebook 1 | 3-4日 | Clarification実装 |
| 3. Notebook 2 | 3-4日 | Setup基本実装 |
| 4. Notebook 3 | 4-5日 | Setup Coordinator |
| 5. Notebook 4 | 2-3日 | Validation実装 |
| 6. Notebook 5 | 3-4日 | 全統合 |
| **合計** | **4週間** | **MVP完成** |

---

### プロジェクト特性（deep_researchパターン採用による強化）

#### 開発効率
- ✅ **インタラクティブ開発**: Notebookで即座にテスト
- ✅ **ドキュメント統合**: 説明とコードが一体
- ✅ **再現性**: Notebook実行でsrc/生成
- ✅ **学習コスト低減**: deep_researchの実装パターンを踏襲

#### 技術的特徴
- ✅ **非競合**: CHARMM-GUIと棲み分け（Amber特化）
- ✅ **将来性**: MCP標準でツール永続化、LLM/実行基盤の更新に強い
- ✅ **拡張性**: LangGraphのモジュラー設計
- ✅ **標準準拠**: LangChain 1.0 + LangGraph 1.0の最新パターン

#### 品質保証
- ✅ **Structured Output**: Pydanticスキーマで決定を明示化
- ✅ **決定ログ**: 全てのツール選択理由を記録
- ✅ **チェックポイント**: 中断・再開・time-travel可能
- ✅ **QC統合**: 物理化学的一貫性チェック

---

### 参考資料

#### 必読
1. **deep_research_from_scratch** (このプロジェクトの参考実装)
   - README.md: 全体アーキテクチャ
   - CLAUDE.md: 開発ワークフロー
   - notebooks/: 5つの実装例

2. **LangGraph公式ドキュメント**
   - Command API: https://langchain-ai.github.io/langgraph/concepts/low_level/#command
   - Subgraphs: https://langchain-ai.github.io/langgraph/how-tos/subgraph/
   - Persistence: https://langchain-ai.github.io/langgraph/concepts/persistence/

3. **MCP Protocol**
   - langchain-mcp-adapters: https://github.com/langchain-ai/langchain-mcp
   - MCP Specification: https://modelcontextprotocol.io/
