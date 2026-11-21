"""SQLAlchemy persistence models for DockLab jobs and scientific outputs."""
from datetime import datetime, timezone
from sqlalchemy import DateTime, Float, ForeignKey, Integer, String, Text, JSON
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


class Job(Base):
    __tablename__ = "jobs"
    id: Mapped[str] = mapped_column(String(36), primary_key=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))
    status: Mapped[str] = mapped_column(String(24), default="created", index=True)
    protein_path: Mapped[str | None] = mapped_column(Text, nullable=True)
    box: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    settings: Mapped[dict] = mapped_column(JSON, default=dict)
    ligands: Mapped[list] = mapped_column(JSON, default=list)
    qc_report: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    results: Mapped[list] = mapped_column(JSON, default=list)
    progress_total: Mapped[int] = mapped_column(Integer, default=0)
    progress_done: Mapped[int] = mapped_column(Integer, default=0)
    error_message: Mapped[str | None] = mapped_column(Text, nullable=True)
    expires_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
