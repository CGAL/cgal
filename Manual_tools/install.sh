#!/bin/sh

  . ./install.config

install -d $LATEX_CONV_BIN $LATEX_CONV_CONFIG $STYLE_FILES $STYLE_FILES/eps_tabs

cd src
#make clean
make || exit 1
make install
cd ..

cp scripts/index_fix scripts/cc_make_ref_pages scripts/cc_ref_wizard scripts/tex2doxy $LATEX_CONV_BIN
cp sty/*.sty $STYLE_FILES
cp sty/eps_tabs/*.pdf $STYLE_FILES/eps_tabs
# cp  sty/eps_tabs_grey/*.eps $STYLE_FILES/eps_tabs_grey # sty/eps_tabs_grey/*.pdf

echo ""
echo "================================================"
echo "Manual_tools successfully installed. Do not forget to update \$TEXINPUTS:"
echo "export TEXINPUTS=\".:$STYLE_FILES:\$TEXINPUTS\""
echo "================================================"
echo "see also src/INSTALLATION for further necessary environment variables, e.g.: "
echo "export LATEX_CONV_CONFIG=\".:$LATEX_CONV_CONFIG:\$LATEX_CONV_CONFIG\"" 
echo "export LATEX_CONV_BIN=\".:$LATEX_CONV_BIN:\$LATEX_CONV_BIN\"" 
echo "================================================"

