texwfiles=$(wildcard *.texw)
pyfiles=$(patsubst %.texw,%.py,$(texwfiles))
pdffiles=$(patsubst %.texw,%.pdf,$(texwfiles))

all: pdfs pys

pys: $(pyfiles)

pdfs: $(pdffiles)

%.py: %.texw
	ptangle $<

%.tex: %.texw
	pweave -f texminted $<

%.pdf: %.tex
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

M2.pdf: bapolar.py
M3.pdf: bapolar.py
det3.pdf: bapolar.py
