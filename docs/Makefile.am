# Documentation for dtq.
#
# The pre-built qd.pdf and the quick-reference guides are installed
# under $(docdir).  The LaTeX sources and figures are distributed so
# the paper can be rebuilt with "make doc" (requires pdflatex and
# bibtex); a normal "make" never invokes LaTeX.

dist_doc_DATA = qd.pdf DTQ_QUICKREF.en.md DTQ_QUICKREF.ja.md

EXTRA_DIST = \
  qd.tex qd.bib \
  nine-two-sum.eps nine-two-sum.fig \
  normal_sum_prod.eps normal_sum_prod.fig \
  qd_add.eps qd_add.fig \
  qd_add_proof.eps qd_add_proof.fig \
  qd_add_qd_d.eps qd_add_qd_d.fig \
  qd_mul_accum.eps qd_mul_accum.fig \
  qd_mul_qd_d.eps qd_mul_qd_d.fig \
  quick-two-sum.eps quick-two-sum.fig \
  six-three-sum.eps six-three-sum.fig \
  three-sum-2.eps three-sum-2.fig \
  three-sum-3.eps three-sum-3.fig \
  three-sum.eps three-sum.fig \
  two-prod.eps two-prod.fig \
  two-sum.eps two-sum.fig

doc:
	cd $(srcdir) && \
	pdflatex qd.tex && \
	bibtex qd && \
	pdflatex qd.tex && \
	pdflatex qd.tex

.PHONY: doc
