ALNS=dnmt1-hym.clu dnmt2-hym.clu dnmt3-hym.clu mbd-hym.clu tet2-hym.clu \
     dnmt1-vert.clu dnmt2-vert.clu dnmt3-vert.clu mbd-vert.clu

all:		$(ALNS)

clean:		
		rm -f $(ALNS)

%.clu:		%.faa
		clustalo --in $< --outfmt clu --wrap 100 > $@
