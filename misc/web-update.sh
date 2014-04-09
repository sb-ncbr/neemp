#!/bin/bash

original_name="neemp-latest.tar.gz"

last_commit=$(git log -1 --pretty=format:"%ci: %s")

user="xracek"
server="aisa.fi.muni.cz"
dir="/home/xracek/public_html/neemp"
webdir="http://fi.muni.cz/~$user/neemp"

cd .. && git archive --worktree-attributes --format=tar.gz HEAD > $original_name

scp "$original_name" $user@$server:$dir
ssh $user@$server "sed -i \"s/^Latest change: .*/Latest change: $last_commit/\" $dir/index.html"
