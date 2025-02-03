let
  pkgs = import <nixpkgs> {};
  pkgs_list = with pkgs.python312.pkgs; 
  	      [
	      pip
	      numpy
	      matplotlib
	      networkx];
  libsbmlCustom = pkgs.libsbml.overrideAttrs (oldAttrs: {
    cmakeFlags = (oldAttrs.cmakeFlags or []) ++ [
		      #"-DWITH_GROUPS=ON"
		      #"-DWITH_COMP=ON"
		      #"-DWITH_FBC=ON"
		      #"-DWITH_LAYOUT=ON"
			"-DWITH_STABLE_PACKAGES=ON"
		    ];
     #withPython = true;
  });
  #mylibsbml = pkgs.callPackage ./libsbml/libsbml.nix {};
in pkgs.mkShell {
  buildInputs = [
    pkgs.python3
    pkgs.libxml2
    pkgs.libz
    #mylibsbml
  ] ++ pkgs_list;
  shellHook = ''
    # Tells pip to put packages into $PIP_PREFIX instead of the usual locations.
    # See https://pip.pypa.io/en/stable/user_guide/#environment-variables.
    export PIP_PREFIX=$(pwd)/_build/pip_packages
    export PYTHONPATH="$PIP_PREFIX/${pkgs.python3.sitePackages}:$PYTHONPATH"
    export PATH="$PIP_PREFIX/bin:$PATH"
    unset SOURCE_DATE_EPOCH
  '';
}
