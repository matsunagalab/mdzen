"""
LM Studio integration for local LLM execution.

Provides OpenAI-compatible API client for LM Studio local server.
"""

import os
import logging
from typing import Optional
from openai import OpenAI

logger = logging.getLogger(__name__)


class LMStudioClient:
    """LM Studio OpenAI互換APIクライアント
    
    LM Studioで実行中のローカルLLMと通信するためのクライアント。
    プライバシー保護のため、データは全てローカルに留まる。
    
    環境変数:
        LM_STUDIO_BASE_URL: LM StudioサーバーURL (デフォルト: http://localhost:1234/v1)
        LM_STUDIO_MODEL: 使用するモデル名 (デフォルト: gpt-oss-20b)
    
    Example:
        >>> client = LMStudioClient()
        >>> if client.is_available():
        ...     response = await client.complete("What is molecular dynamics?")
        ...     print(response)
    """
    
    def __init__(
        self,
        base_url: Optional[str] = None,
        model: Optional[str] = None
    ):
        """Initialize LM Studio client
        
        Args:
            base_url: LM StudioサーバーURL（Noneの場合は環境変数から取得）
            model: 使用するモデル名（Noneの場合は環境変数から取得）
        """
        self.base_url = base_url or os.getenv(
            "LM_STUDIO_BASE_URL",
            "http://localhost:1234/v1"
        )
        self.model = model or os.getenv(
            "LM_STUDIO_MODEL",
            "gpt-oss-20b"
        )
        
        # OpenAI clientの初期化（LM Studioは認証不要なのでダミーキー）
        self.client = OpenAI(
            base_url=self.base_url,
            api_key="lm-studio"
        )
        
        logger.info(f"LM Studio client initialized: {self.base_url}, model={self.model}")
    
    async def complete(
        self,
        prompt: str,
        system: Optional[str] = None,
        temperature: float = 0.3,
        max_tokens: int = 2000
    ) -> str:
        """LLM completion実行
        
        Args:
            prompt: ユーザープロンプト
            system: システムプロンプト（オプション）
            temperature: サンプリング温度（0.0-1.0、低いほど決定的）
            max_tokens: 最大トークン数
        
        Returns:
            LLMの応答テキスト
        
        Raises:
            Exception: LM Studioサーバーとの通信エラー
        """
        messages = []
        
        if system:
            messages.append({"role": "system", "content": system})
        
        messages.append({"role": "user", "content": prompt})
        
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                temperature=temperature,
                max_tokens=max_tokens
            )
            
            content = response.choices[0].message.content
            logger.debug(f"LLM response: {content[:100]}...")
            return content
            
        except Exception as e:
            logger.error(f"LM Studio API error: {e}")
            raise
    
    def complete_sync(
        self,
        prompt: str,
        system: Optional[str] = None,
        temperature: float = 0.3,
        max_tokens: int = 2000
    ) -> str:
        """LLM completion実行（同期版）
        
        Args:
            prompt: ユーザープロンプト
            system: システムプロンプト（オプション）
            temperature: サンプリング温度
            max_tokens: 最大トークン数
        
        Returns:
            LLMの応答テキスト
        """
        messages = []
        
        if system:
            messages.append({"role": "system", "content": system})
        
        messages.append({"role": "user", "content": prompt})
        
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                temperature=temperature,
                max_tokens=max_tokens
            )
            
            return response.choices[0].message.content
            
        except Exception as e:
            logger.error(f"LM Studio API error: {e}")
            raise
    
    def is_available(self) -> bool:
        """LM Studioサーバー接続確認
        
        Returns:
            True if LM Studio server is running and accessible
        """
        try:
            # モデルリストを取得して接続確認
            models = self.client.models.list()
            logger.info(f"LM Studio server available. Models: {len(models.data)}")
            return True
        except Exception as e:
            logger.warning(f"LM Studio server not available: {e}")
            return False
    
    def get_available_models(self) -> list[str]:
        """利用可能なモデルのリストを取得
        
        Returns:
            モデル名のリスト
        """
        try:
            models = self.client.models.list()
            return [model.id for model in models.data]
        except Exception as e:
            logger.error(f"Failed to get models: {e}")
            return []


# シングルトンインスタンス（簡単に使えるように）
_default_client: Optional[LMStudioClient] = None


def get_default_client() -> LMStudioClient:
    """デフォルトのLM Studioクライアントを取得
    
    Returns:
        シングルトンLMStudioClientインスタンス
    """
    global _default_client
    if _default_client is None:
        _default_client = LMStudioClient()
    return _default_client

