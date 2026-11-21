"""Celery application for asynchronous docking work."""
from celery import Celery
from backend.app.config import settings

celery_app = Celery("docklab", broker=settings.CELERY_BROKER_URL, backend=settings.CELERY_RESULT_BACKEND)
celery_app.conf.update(task_track_started=True, task_serializer="json", result_serializer="json", accept_content=["json"], timezone="UTC")
celery_app.autodiscover_tasks(["backend.app.workers"])
