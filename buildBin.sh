#! /bin/sh
arch=`uname -m`
version=1.0.0-1
# Note: run from compilation dir (e.g. build/)
strip grapes/grapes
tar cvzf grapes-${arch}-bin-static-${version}.tar.gz grapes/grapes

