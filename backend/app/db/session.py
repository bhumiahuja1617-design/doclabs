"""Database wiring; SQLite supports a local no-Docker development run."""
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from backend.app.config import settings
from .models import Base

url = settings.DATABASE_URL.replace("sqlite+aiosqlite://", "sqlite://").replace("postgresql+asyncpg://", "postgresql+psycopg://")
engine = create_engine(url, pool_pre_ping=True)
SessionLocal = sessionmaker(engine, expire_on_commit=False)


def init_db() -> None:
    Base.metadata.create_all(engine)
