# 目录组织

- **Solver/**
  旅行货郎担问题的 Visual C++ 项目.
  核心算法.

- **Protocol/**
  - TSP.proto
    旅行货郎担问题的输入输出数据定义文件.

- **Simulator/**
  提供算法交互接口的 Visual C++ 项目.
  实现批量测试与结果展示等功能.

- **Analyzer/**
  旅行货郎担问题算例/结果/日志分析程序的 Visual C# 项目.

- **Deploy/**
  开发阶段程序部署目录, 包含程序运行所需要的所有可执行文件和数据.
  - **Instance/**
    输入数据.
  - **Solution/**
    求解结果.
  - **Visualization/**
    结果可视化.

- **Lib/**
  第三方库.
  - **protobuf**
    由 Google 开发的数据交换格式.
    访问 https://github.com/google/protobuf/releases 下载最新版本.
    - **bin/**
      生成对消息进行序列化与反序列化的代码的工具.
    - **include/**
      头文件. 应被设置为 "附加包含文件".
    - **lib/**
      静态链接库. 应被设置为 "附加库目录".
      该文件夹下的库文件应被设置为 "附加依赖项".



# 编译链接


## 编译链接

1. 自己编译 protobuf 并解压至 `Lib` 目录.
2. 打开根目录下的 Visual Studio 2017 解决方案文件 TSP.sln`, 设置输出选项为 Release x64（这里只配置了Release编译环境）, 生成解决方案.
