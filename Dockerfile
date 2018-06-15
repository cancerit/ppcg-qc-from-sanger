FROM python:3.6-alpine

MAINTAINER yx2@sanger.ac.uk

LABEL uk.ac.sanger.cgp="Cancer Genome Project, Wellcome Trust Sanger Institute" \
      description="tool to extract PPCG defined QC metrics for dockstore.org"

# Add service user
ENV USER=cgp
ENV HOME=/opt/cgp
RUN mkdir -p $HOME
RUN addgroup $USER && \
    adduser -G $USER -D -h $HOME $USER
WORKDIR $HOME

# Temporarily install alpine linux dependencies for building python, ruby and alpine linux dependencies
RUN apk add --no-cache --virtual build-deps build-base linux-headers
# Install alpine linux dependencies
RUN apk add --no-cache bash ruby

# Removes the need for rdoc when installing mdl
RUN echo 'gem: --no-document' > /etc/gemrc
# Install ruby requirement for parsing markdown files
RUN gem install mdl

WORKDIR $HOME
COPY requirements.txt ./
RUN pip install -r requirements.txt

# Clear out alpine linux build dependencies
RUN apk del build-deps

# Copy over source code
COPY . .

# Build core tools
RUN python setup.py develop

RUN chown -R $USER:$USER .

# Become the final user
USER cgp

CMD ls
