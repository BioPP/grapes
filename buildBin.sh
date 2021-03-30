#! /bin/sh
arch=`uname -m`
version=1.1.1-1
# Note: run from compilation dir (e.g. build/)
strip grapes/grapes
strip grapes/multi_grapes
tar cvzf grapes-${arch}-bin-static-${version}.tar.gz grapes/grapes grapes/multi_grapes

