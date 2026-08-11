ACLOCAL_AMFLAGS = -I m4
SUBDIRS = config src include tests docs

dist_doc_DATA = README README.md

EXTRA_DIST = LICENSE

doc:
	$(MAKE) -C docs doc

changelog:
	@if test -d "$(top_srcdir)/.git"; then \
	  (cd $(top_srcdir) && git log) >ChangeLog; \
	else \
	  echo "Not a git repository; skipping changelog generation."; \
	fi

dist-hook:
	@if test -d "$(top_srcdir)/.git"; then \
	  rm -f $(distdir)/ChangeLog && \
	  (cd $(top_srcdir) && git log) >$(distdir)/ChangeLog; \
	else \
	  echo "Not a git repository; skipping ChangeLog regeneration in $(distdir)."; \
	  if test -f $(top_srcdir)/ChangeLog && test ! -f $(distdir)/ChangeLog; then \
	    cp $(top_srcdir)/ChangeLog $(distdir)/ChangeLog; \
	  fi; \
	fi

cpp-demo:
	$(MAKE) -C tests demo

demo: cpp-demo

time:
	$(MAKE) -C tests time

bin_SCRIPTS=qd-config

.PHONY: changelog cpp-demo demo time doc
