# Use an official Fedora as a parent image for the build stage
FROM fedora:latest AS sources_deps

# Set environment variables to non-interactive
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN dnf update -y && dnf install -y \
    wget \
    make \
    gcc \
    gcc-c++ \
    patch \
    cmake \
    bison \
    flex \
    unzip \
    python3 \
    && dnf clean all

# Copy the patch file to the build context
COPY cgal-NO_ADDITIONAL_DETAILS.patch .

FROM sources_deps AS build

# Build and install Doxygen from sources
ARG DOXYGEN_VERSION=1.9.6
ARG MAKEFLAGS=-j$(nproc)
RUN if [ -n "$DEBUG"];then set -x && make --version && ls -lZ /tmp && id; fi \
    && DOXYGEN_VERSION_UNDERSCORE=$(echo ${DOXYGEN_VERSION} | sed 's/\./_/g') \
    && wget https://github.com/doxygen/doxygen/archive/refs/tags/Release_${DOXYGEN_VERSION_UNDERSCORE}.zip \
    && unzip Release_${DOXYGEN_VERSION_UNDERSCORE}.zip \
    && cd doxygen-Release_${DOXYGEN_VERSION_UNDERSCORE} \
    && patch -p1 < ../cgal-NO_ADDITIONAL_DETAILS.patch \
    && mkdir build \
    && cd build \
    && cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release .. \
    && cmake --build . \
    && cmake --install . \
    && mkdir -p /usr/local/share/doc/doxygen && cp ../LICENSE /usr/local/share/doc/doxygen/LICENSE.TXT \
    && cd ../.. \
    && rm -rf doxygen-Release_${DOXYGEN_VERSION_UNDERSCORE} Release_${DOXYGEN_VERSION_UNDERSCORE}.zip

# Use a smaller base image for the final stage
FROM fedora:latest

# Install necessary runtime dependencies
RUN set -x \
    && dnf update -y && dnf install -y graphviz 'perl(Getopt::Std)' tex-bibtex cmake python3-lxml python3-pyquery \
    && dnf clean all

# Copy Doxygen from the build stage
COPY --from=build /usr/local/bin/doxygen /usr/local/bin
COPY --from=build /usr/local/share/doc/doxygen/LICENSE.TXT /usr/local/share/doc/doxygen/LICENSE.TXT
RUN doxygen --version
