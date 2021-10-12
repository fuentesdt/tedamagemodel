.PHONY: doc rtf
PAPER=damagemodel
all: 
	pdflatex $(PAPER).tex
bib:
	bibtex  $(PAPER)
clean:
	rm  $(PAPER).aux  $(PAPER).bbl  $(PAPER).blg  $(PAPER).log  $(PAPER).out  
