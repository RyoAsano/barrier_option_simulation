export barrier_project_root_dir=~/my_repo/barrier_option_simulation/time_changed_barrier_option2
sed -i "" -e 's/\-isysroot \/Library\/Developer\/CommandLineTools\/SDKs\/MacOSX10.14\.sdk//g' $barrier_project_root_dir/build/compile_commands.json
