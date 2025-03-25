# Use a base image for Intel architecture (amd64)
FROM --platform=linux/amd64 ubuntu:20.04

# Update package list and install build-essential, which includes make, gcc, etc.
RUN apt-get update && apt-get install -y build-essential make

# Set the working directory for the Makefile
WORKDIR /app

# Copy your source code and Makefile into the container
COPY . /app

# Run `make` command (you can also specify target if needed)
RUN make