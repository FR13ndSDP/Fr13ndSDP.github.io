# Happy FOAMing on VScode




![badge](https://img.shields.io/badge/Modern-Editor-green)

使用现代化的编辑器vscode进行OpenFOAM编程的配置，大大提升编程体验。
<!--more-->

### 说明
- 首先搜索扩展`c/c++`安装微软的的C++扩展

- 打开工作文件夹，即你的OpenFOAM求解器存放的文件夹，按住`shift+ctrl+p`开启命令框，输入`c/c++:conf`后选择如图的第二条命令，生成`.vscode/c_cpp_properties.json`

{{< figure src="/images/vscode-openfoam/step1.png" >}}

- 打开`Make/options`，可以看到所有需要包含的头文件，以`-I`开头

{{< figure src="/images/vscode-openfoam/step2.png" >}}

值得注意的是，其中的`$(LIB_SRC)`所指的地址就是OpenFOAM源代码路径。

如果你不知道这个路径，终端下执行

```
src
pwd
```
打印出的内容就是`$LIB_SRC`

**值得注意的是，`Make/options`中的头文件并不完整，还需要加入以下两个头文件：**
```
/opt/OpenFOAM/OpenFOAM-v1912/USERNAME/src/OpenFOAM/lnInclude
/opt/OpenFOAM/OpenFOAM-v1912/USERNAME/src/OSspecific/POSIX/lnInclude
```
请自行搜索以上两个头文件在你的设备中的位置。

- 接着把完整路径放到C++配置文件的`includePath`中：

{{< figure src="/images/vscode-openfoam/step3.png" >}}

注意每条路径用`" "`括起来，行末加`,`（最后一行路径除外）

- 接下来还有一些小工作要做，先`wmake`求解器，你会看到下面的输出：

{{< figure src="/images/vscode-openfoam/step4.png" >}}

第一个参数`g++`说明`wmake`使用g++编译器，因此将g++的路径替换`compilerPath`中的内容，不知道g++路径？终端里`which g++`得到。后面`-std = c++11 -m64 ... -ftemplate-depth-100`是`compilerArgs`，
大部分是规定了编译过程错误信息输出以及优化信息，其中，`-W`开头的参数规定了编译过程中错误信息如何输出，对于我们使用编辑器没有作用，`-std = c++ 11 -m64 -O3`这几个参数则规定了编译使用的C++标准和优化信息，可以保留，而
以`-D`参数开头的是一些宏定义，他们有的规定了所使用的精度等级，这些定义十分重要，如果不加会引起类型错误。最后`-ftemplate-depth-100`规定了实例化搜寻深度，默认值为900，所以如果不加应该也没有问题，但是对性能的影响尚不清楚。
如果你觉得太麻烦，你可以将他们全都加到`compilerArgs`里

{{< figure src="/images/vscode-openfoam/step5.png" >}}

- 最终效果：

{{< figure src="/images/vscode-openfoam/step6.png" >}}

自动补全，代码提示、跳转都正常

-----

### NOTE
- 由于OpenFOAM”特殊“的编程风格，常在头文件中写代码段而不是类声明等常规写法，会导致报错，比如头文件中`runTime`的`undefined`错误提示和`Info`的`ambigous`提示，目前没有解决方法。如果你不想看到这些错误提示，可以设置忽略，具体方法在网上很容易找到

- 最近有个叫Volker Weissmann 的哥们在github上发起一个pull request, 介绍了一个叫做`ccls`的东西，似乎能够解决上面的问题，不过我还没尝试。[pull request](https://github.com/OpenFOAM/OpenFOAM-dev/pull/31)

  {{< figure src="/images/vscode-openfoam/mail.png" >}}

- OpenFOAM的`.C .H`后缀可能使得Vscode识别语言类型错误，可以在`.vscode`下创建`settings.json`设置类型绑定，具体方法在网上很容易找到

- 由于C++扩展进行intellisense会占用较大的存储空间，每打开一个工程，就会相应生成一个100Mb左右的预编译文件，但是对于intellisense没有作用。可以定期清理`~/.cache/vscode-cpptools/ipch`下的缓存文件来释放空间，或者在设置中把`intellisense cache size`设为0禁用。

- 如果是wsl用户，可以使用wsl插件在vscode中调用bash，*上述操作均在vscode连接Ubuntu后进行*

- 如果遇到问题请[点我][1]发起`issue`

----

### LICENSE

本内容基于CC-BY-4.0 协议，转载请注明出处

[1]:https://gitee.com/fr13ndsdp/Vscode-OpenFOAM/issues/new
