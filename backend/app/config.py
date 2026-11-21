from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field

class Settings(BaseSettings):
    PROJECT_NAME: str = "DockLab"
    ENV: str = "development"
    DEBUG: bool = True
    
    # Paths
    BASE_DIR: Path = Path(__file__).resolve().parent.parent.parent
    DATA_DIR: Path = Path(__file__).resolve().parent.parent.parent / "data"
    
    # Database
    DATABASE_URL: str = "sqlite+aiosqlite:///./docklab.db"
    TEST_DATABASE_URL: str = "sqlite+aiosqlite:///./test.db"
    
    # Redis & Queue
    REDIS_URL: str = "redis://redis:6379/0"
    CELERY_BROKER_URL: str = "redis://redis:6379/0"
    CELERY_RESULT_BACKEND: str = "redis://redis:6379/1"
    
    # Docking & Vina
    DEFAULT_EXHAUSTIVENESS: int = 8
    MAX_EXHAUSTIVENESS: int = 32
    DEFAULT_POSES: int = 9
    DEFAULT_ENERGY_RANGE: float = 3.0
    DEFAULT_VINA_TIMEOUT_SECONDS: int = 600
    
    # Reliability Settings
    RELIABILITY_DEFAULT_TOP_K: int = 10
    RELIABILITY_MAX_TOP_K: int = 20
    RELIABILITY_DEFAULT_SEEDS: int = 5
    RELIABILITY_MAX_SEEDS: int = 10
    RELIABILITY_HIGH_CONSISTENCY: float = 0.80
    RELIABILITY_HIGH_MAX_SD: float = 0.30
    RELIABILITY_MODERATE_CONSISTENCY: float = 0.50
    RELIABILITY_MODERATE_MAX_SD: float = 0.70
    
    # Quotas & Limits
    MAX_UPLOAD_SIZE_MB: int = 10
    MAX_MOLECULES_PER_JOB: int = 200
    MAX_BOX_SIZE_ANGSTROMS: float = 30.0
    WARN_BOX_SIZE_ANGSTROMS: float = 25.0
    JOB_RETENTION_DAYS: int = 7
    
    # CORS
    CORS_ORIGINS: list[str] = [
        "http://localhost:5173",
        "http://localhost:3000",
        "http://127.0.0.1:5173",
    ]
    
    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore"
    )

settings = Settings()

