SHELL := /bin/bash
DOXYGEN_VERSIONS := 1.12.0 1.11.0 1.10.0 1.9.8 1.9.7 1.9.6

.PHONY: all build-% push-% build push

all: build
	@echo "Use `$(MAKE) push` to push the images to the registry."

build-%:
	@echo "MAKEFLAGS: $(MAKEFLAGS)"
	@echo "Building Doxygen version $*..."
	if [ "$$(getenforce || true)" == "Enforcing" ]; then Z=:z; else Z=; fi; \
	F="$(MAKEFLAGS)"; F=$${F##*fifo:}; F=$${F%% *}; \
	if [ -p "$$F" ]; then echo "The GNU make FIFO file exists:"; ls -l $$F; VOLUME_ARGS="-v $$F:$$F$$Z"; echo -- $$VOLUME_ARGS; fi; \
	podman build --build-arg DOXYGEN_VERSION=$* --build-arg "MAKEFLAGS=$(MAKEFLAGS)" $$VOLUME_ARGS -t cgal/doxygen:$* .

push-%: build-%
	@echo "Pushing cgal/doxygen:$*..."
	podman push cgal/doxygen:$*

build: $(foreach version,$(DOXYGEN_VERSIONS),build-$(version))

push: $(foreach version,$(DOXYGEN_VERSIONS),push-$(version))
