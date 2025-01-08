说明：

​	bin目录中有一个名为“icts”的可执行文件，如果需要重新编译，则在根目录下使用 make 指令即可，编译成功会在 bin 目录中生成 icts 可执行文件

执行方式：

​	./icts -def_file ${电路文件全路径} -constraint ${约束文件全路径} -output ${solution文件全路径}

示例：

​	./icts -def_file /home/public/case1/problem.def -constraint /home/public/case1/constraints.txt -output /home/public/case1/solution.def

