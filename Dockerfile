FROM python:3.6

USER root

WORKDIR /app

ADD . /app

COPY cleanup-cron /etc/cron.d/cleanup-cron

EXPOSE 80

RUN pip install -r requirements.txt \
&& apt-get update && apt-get install --yes cron \
&& chmod 0644 /etc/cron.d/cleanup-cron \
&& chmod 0744 cleanup.sh \
&& crontab /etc/cron.d/cleanup-cron

CMD cron && gunicorn -w 4 -b :80 --preload lavaruins:server
