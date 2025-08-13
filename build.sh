#!/Bin/bash

#
mkdir -p ~/Documents/Centerpoint/build
pushd  ~/Documents/Centerpoint/build

#oprimized build:

#g++  ../smalltest.cpp  -o centerpoint -DPOLYMAKE_VERSION=413 -fPIC -pipe -std=c++14 -Wno-logical-op-parentheses -Wno-shift-op-parentheses -Wno-mismatched-tags -Wno-unused-local-typedef -Wno-error=unneeded-internal-declaration -Wshadow -Wconversion -Wno-sign-conversion -Wzero-as-null-pointer-constant -DPOLYMAKE_WITH_FLINT -O3 -DPOLYMAKE_DEBUG=0 -lpolymake -lpolymake-apps -lc++ -lflint -lmpfr -lgmp -lpthread -I/opt/homebrew/Cellar/polymake/4.13_3/include -I/opt/homebrew/Cellar/polymake/4.13_3/include/polymake/external -I/opt/homebrew/include -L/opt/homebrew/Cellar/polymake/4.13_3/lib -L/opt/homebrew/lib

#debug build:
g++  ../smalltest.cpp  -o centerpoint -DPOLYMAKE_VERSION=413 -fPIC -pipe -std=c++14 -Wno-logical-op-parentheses -Wno-shift-op-parentheses -Wno-mismatched-tags -Wno-unused-local-typedef -Wno-error=unneeded-internal-declaration -Wshadow -Wconversion -Wno-sign-conversion -Wzero-as-null-pointer-constant -Wno-zero-as-null-pointer-constant -Wno-null-dereference -DPOLYMAKE_WITH_FLINT -g -DPOLYMAKE_DEBUG=1 -lpolymake -lpolymake-apps -lc++ -lflint -lmpfr -lgmp -lpthread -lgurobi_c++ -lgurobi120 -I/opt/homebrew/Cellar/polymake/4.13_3/include -I/opt/homebrew/Cellar/polymake/4.13_3/include/polymake/external -I/opt/homebrew/include -I/Library/gurobi1201/macos_universal2/include -L/opt/homebrew/Cellar/polymake/4.13_3/lib -L/opt/homebrew/lib -L/Library/gurobi1201/macos_universal2/lib

#clang -g -fuse-ld=lld ../sim86.cpp -o sim86_clang_debug.exe
# do -g for debug build
# -L/Library/gurobi1201/macos_universal2/lib
# -I/Library/gurobi1201/macos_universal2/include

#For using the Assert makro: -Wno-zero-as-null-pointer-constant -Wno-null-dereference

#clang++ -std=c++14 -O3 -fPIC -pipe -DPOLYMAKE_VERSION=413 -DPOLYMAKE_DEBUG=0 -DPOLYMAKE_WITH_FLINT -Wno-missing-template-arg-list-after-template-kw -Wno-logical-op-parentheses -Wno-shift-op-parentheses -Wno-mismatched-tags -Wno-unused-local-typedef -Wno-error=unneeded-internal-declaration -Wshadow -Wconversion -Wno-sign-conversion -Wzero-as-null-pointer-constant -I/opt/homebrew/Cellar/polymake/4.13_3/include -I/opt/homebrew/Cellar/polymake/4.13_3/include/polymake/external -I/opt/homebrew/include -L/opt/homebrew/Cellar/polymake/4.13_3/lib -L/opt/homebrew/lib ../smalltest.cpp -o smalltest -lpolymake -lpolymake-apps -lflint -lmpfr -lgmp -lc++ -lpthread
popd
