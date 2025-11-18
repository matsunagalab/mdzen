import sys
import asyncio
from pathlib import Path
from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain.chat_models import init_chat_model
from langchain.agents import create_agent

async def main():
    ### MCP tool setup
    project_root = Path(__file__).parent.parent
    server_config = {
        "structure": {
            "transport": "stdio",
            "command": sys.executable,
            "args": [str(project_root / "servers" / "structure_server.py")],
        }
    }

    ### create_agent
    mcp_md_client = MultiServerMCPClient(server_config)
    tools = await mcp_md_client.get_tools()

    model = init_chat_model("ollama:gpt-oss:20b")

    md_agent = create_agent(
        model,
        tools, 
        system_prompt="""You are a smart AI-assistant. """
    )

    ### chat with agent
    print("🧪 Testing agent response...")

    result = await md_agent.ainvoke({
        "messages": [{"role": "user", "content": "人間とはなんでしょう？1回目回答した後で、自分の回答を振り返り、よりわかりやすくなるように改善すべき点は改善して、2回目の回答をしてください。それを5回目まで繰り返してください"}]
    })

    print("\n🤖 Agent says:")
    for msg in result["messages"]:
        if msg.type == "ai" and msg.content:
            print(msg.content)

    print("\n✅ Test done!")

if __name__ == "__main__":
    asyncio.run(main())