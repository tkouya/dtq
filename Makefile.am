ACLOCAL_AMFLAGS = -I m4
SUBDIRS = config src include tests

BUILT_SOURCES = docs/qd.pdf
dist_doc_DATA = README README.md docs/qd.pdf \
                docs/DTQ_QUICKREF.en.md \
                docs/DTQ_QUICKREF.ja.md
dist_noinst_DATA = docs/Makefile \
  docs/qd.tex docs/qd.bib \
  docs/nine-two-sum.eps docs/nine-two-sum.fig \
  docs/normal_sum_prod.eps docs/normal_sum_prod.fig \
  docs/qd_add.eps docs/qd_add.fig \
  docs/qd_add_proof.eps docs/qd_add_proof.fig \
  docs/qd_add_qd_d.eps docs/qd_add_qd_d.fig \
  docs/qd_mul_accum.eps docs/qd_mul_accum.fig \
  docs/qd_mul_qd_d.eps docs/qd_mul_qd_d.fig \
  docs/quick-two-sum.eps docs/quick-two-sum.fig \
  docs/six-three-sum.eps docs/six-three-sum.fig \
  docs/three-sum-2.eps docs/three-sum-2.fig \
  docs/three-sum-3.eps docs/three-sum-3.fig \
  docs/three-sum.eps docs/three-sum.fig \
  docs/two-prod.eps docs/two-prod.fig \
  docs/two-sum.eps docs/two-sum.fig

doc: docs/qd.pdf

docs/qd.pdf:
	$(MAKE) -C docsqd.pdf

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
