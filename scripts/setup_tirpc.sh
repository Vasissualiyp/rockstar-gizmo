#!/usr/bin/env bash

tirpc_url_devel="https://buildlogs.centos.org/centos/7/os/x86_64-latest/Packages/libtirpc-devel-0.2.4-0.3.el7.x86_64.rpm"
tirpc_fname_devel="libtirpc-devel-0.2.4-0.3.el7.x86_64.rpm"
tirpc_dir_devel="libtirpc-devel"

tirpc_url="https://buildlogs.centos.org/centos/7/os/x86_64-latest/Packages/libtirpc-0.2.4-0.3.el7.x86_64.rpm"
tirpc_fname="libtirpc-0.2.4-0.3.el7.x86_64.rpm"
tirpc_dir="libtirpc"

# Setup devel version
mkdir $tirpc_dir
cd $tirpc_dir
wget $tirpc_url
rpm2cpio $tirpc_fname | cpio -idmv
cd ..

# Setup regular version
mkdir $tirpc_dir_devel
cd $tirpc_dir_devel
wget $tirpc_url_devel
rpm2cpio $tirpc_fname_devel | cpio -idmv
