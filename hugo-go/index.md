# HUGO搭建博客


使用HUGO+Github Pages搭建个人博客

<!--more-->

## 安装HUGO

由于WSL2的LocalHost和主机不一样，导致没办法本地预览，我放弃了在WSL上安装hugo，而是安装在Windows上，然后通过WSL调用`hugo.exe`来完成工作流。

在 [HUGO的github仓库](https://github.com/gohugoio/hugo/) 下载windows Release版本，由于下载的很慢，我使用github镜像 [link](https://github.wuyanzheshui.workers.dev/) 来下载然后把hugo.exe的路径添加到系统环境变量。终端里执行

```bash
hugo.exe version
```

打印版本信息，一切正常。

## 配置
```bash
hugo new site /path/to/site
```

后，将建立一个文件夹，结构如下

```
.
├── archetypes
├── config.toml
├── content
├── data
├── layouts
├── resources
├── static
└── themes
```

首先要安装一个主题，我选择了`LoveIt`主题，[仓库地址](https://github.com/dillonzq/LoveIt/ )。

```bash
git clone git@github.com:dillonzq/LoveIt.git themes/LoveIt
```

然后需要编辑`config.toml`，我直接把`/themes/LoveIt/exampleSite`里的给复制了过来，然后略加一点修改就行了。注释非常详尽，参照着改就行了。[link](https://hugoloveit.com/)的文档就更加详尽了，以后遇到坑再看吧。

`static`用来放图片，`aechetypes`里的`default.md`在每次创建文章的时候自动加一个模板。

## 发布内容

```bash
hugo new posts/post.md
```
确认要发布后，然后把`draft:true` 改为`false`

然后

```bash
hugo server
```

就可以在1313端口预览网页了。

## 部署到Github Pages

1. 创建github仓库，名字为`UserName.github.io`
2. 进入`public`文件夹，执行`git init`，`git remote add xxx`
3. 将`public`文件夹内容`push`到github


