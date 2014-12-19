VERSION := 0.14.12.19
TARGET  := HGPP_${VERSION}
BIN_DIR = bin
# Copy script into bin
HGPP:
	echo "Building HGPP"
	mkdir -p $(BIN_DIR)
	cd $(BIN_DIR); tar -zxvf ../src/scripts/scripts.tar.gz
	cd src/hgpp_anchor_fosmid; $(MAKE)
	cd src/hgpp_errcorr; $(MAKE)
	cd src/hgpp_cnn; $(MAKE)
	cd extensions/ASE; $(MAKE)
	cd extensions/ASM/tDMR_detection; $(MAKE)

clean:
	rm -rf $(BIN_DIR)
	cd src/hgpp_anchor_fosmid; $(MAKE) clean
	cd src/hgpp_errcorr; $(MAKE) clean 
	cd src/hgpp_cnn; $(MAKE) clean 
	cd extensions/ASE; $(MAKE) clean
	cd extensions/ASM/tDMR_detection; $(MAKE) clean

.PHONY: clean
