{
  description = "C++ project setup for N-Body problem";

  inputs = {
    # Pointing to the current stable release of nixpkgs. You can
    # customize this to point to an older version or unstable if you
    # like everything shining.
    #
    # E.g.
    #
    # nixpkgs.url = "github:NixOS/nixpkgs/unstable";
    nixpkgs.url = "github:NixOS/nixpkgs/23.11";

    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, ... }@inputs: inputs.utils.lib.eachSystem [
    # Add the system/architecture you would like to support here. Note that not
    # all packages in the official nixpkgs support all platforms.
    "x86_64-linux" "i686-linux" "aarch64-linux" "x86_64-darwin"
  ] (system: let
    pkgs = import nixpkgs {
      inherit system;

      # Add overlays here if you need to override the nixpkgs
      # official packages.
      overlays = [];
      
      # Uncomment this if you need unfree software (e.g. cuda) for
      # your project.
      #
      # config.allowUnfree = true;
    };
  in {
    devShells.default = pkgs.mkShell.override { stdenv = pkgs.gcc13Stdenv; } rec {
      # Update the name to something that suites your project.
      name = "Comp-Astro";

      packages = with pkgs; [
        # Development Tools
		eigen # eigen2?
		mathgl
        cmake
        cmakeCurses
		llvmPackages.openmp
		tbb
		clang-tools_16
		clang 
		bear
		gdb
		# cppcheck
      ];

      # Setting up the environment variables you need during
      # development.
      shellHook = ''
		echo "Welcome to the Dev shell for N-body Problem"
		'';
    };

  });
}
