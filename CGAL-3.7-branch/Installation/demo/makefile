all:
	@for DIR in *; do \
	  if test -d "$$DIR"; then \
	    ( cd "$$DIR" && $(MAKE) ) \
	  fi; \
	done

clean:
	@for DIR in *; do \
	  if test -d "$$DIR"; then \
	    ( cd "$$DIR" && $(MAKE) clean ) \
	  fi; \
	done
