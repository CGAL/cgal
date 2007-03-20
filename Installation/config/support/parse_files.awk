BEGIN { 
    RPKGS=""; PKGS="";
    for (i = 1; i < ARGC; i++) {
	pkg=ARGV[i]; sub(".*-","",pkg);
	RPKGS=pkg " " RPKGS; PKGS=PKGS " " pkg;
    };
    printf "SUPPORT_PKGS='%s'\n",PKGS;
    printf "R_SUPPORT_PKGS='%s'\n",RPKGS;
    FS=" *= *"
}
/(PROVIDES|DESCRIPTION|CXXFLAGS|LDFLAGS|LIBS|REQUIRES|INCOMPATIBLE|STDINCLDIRS|INCLTHING|STDLIBDIRS|LIBTHING|COMPILETESTFLAGS) *=/ {
    sub(".*-","",FILENAME);
    printf "%s_%s='",FILENAME,$1;
    for (i=2; i<NF; i++)
	printf "%s=",$(i);
    printf "%s'\n",$(NF)
}
