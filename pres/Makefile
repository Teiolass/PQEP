all: la_presentazione.pdf

la_presentazione.pdf: pres.tex
	pdflatex pres.tex
	mv pres.pdf la_presentazione.pdf

whole: marticle.cls mstyle.sty rel.tex
	pdflatex rel.tex
	pdflatex rel.tex
	pdflatex rel.tex
	mv rel.pdf la_relazione.pdf
