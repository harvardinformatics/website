FROM python:3.10-slim-bullseye AS pelican

EXPOSE 80
RUN apt update && apt install -y --no-install-recommends git \
  && rm -rf /var/lib/apt/lists/*

ADD requirements.txt requirements.txt
RUN python3 -m pip install --no-cache-dir -r requirements.txt

WORKDIR /app
COPY . .
RUN pelican

FROM httpd:2.4.52-alpine

COPY --from=pelican /app/output/ /usr/local/apache2/htdocs/
