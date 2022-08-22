FROM zldcedars/lizard:latest

COPY helpers.R /etc/R/helpers.R
RUN dir  
ENV PORT=8080
CMD ["/usr/bin/R", "--no-save", "--gui-none", "-f", "/app/run.R"]
