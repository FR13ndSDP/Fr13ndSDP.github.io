# DHT11温湿度传感器——51单片机


DHT11温湿度传感器的使用——C51代码

<!--more-->

### 测量参数

测量范围：20-90% RH 0-50℃

测湿精度：±5% RH

测温精度：±2℃

分辨力：1

### 特性

- 宽电压供电，3V~5.5V
- 数据线接上拉电阻，当长度小于20m，使用5k电阻
- 上电后有1s的不稳定状态
- **DHT11 只有在接收到开始信号后才触发一次温湿度采集**，如果没有接收到主机发送复位信号，DHT11 不主动进行温湿度采集。当数据采集完毕且无开始信号后，DHT11 自动切换到低速模式。

### 时序

{{< figure src="/images/DHT11/sequential.png" title="时序" >}}

#### 1. 起始信号

总线空闲状态为高电平，主机把总线拉低等待DHT11响应；

与MCU相连的SDA数据引脚置为输出模式；

主机把总线拉低至少18毫秒，然后拉高20-40us等待DHT返回响应信号；

#### 2. 读取DHT11响应

SDA数据引脚设为输入模式；

DHT11检测到起始信号后，会将总线拉低80us，然后拉高80us作为响应；

#### 3. 送出数据

DHT11将送出40bit的数据

**格式**: 40bit数据=8位湿度整数+8位湿度小数+8位温度整数+8位温度小数+8位校验，当温湿度数据和等于校验位时校验通过。

似乎小数位全是零，而且整数数据的计算如下例所示

`00011000 00000000` = 24.0

因此通过循环移位接收高八位数据后转为字符串输出即可。

DHT11 在拉高总线 80us 后开始传输数据。每 1bit 数据都以 50us 低电平时隙开始，告诉主机开始传输一位数据了。DHT11 以高电平的长短定义数据位是 0 还是 1，当 50us 低电平时隙过后拉高总线，高电平持续 26~28us 表示数据“0”；持续 70us 表示数据“1”。

当最后 1bit 数据传送完毕后，DHT11 拉低总线 50us，表示数据传输完毕，随后总线由上拉电阻拉高进入空闲状态。

### 代码

```c
#include <reg51.h>
#include <intrins.h>

#define uchar unsigned char
#define uint unsigned int

sbit Data=P1^0;   //定义数据线
uchar rec_dat[12];   //用于显示的接收数据数组

void DHT11_delay_ms(uint z)
{
   uint i,j;
   for(i=z;i>0;i--)
      for(j=110;j>0;j--);
}

void DHT11_start()  // 开始信号
{
   Data=1;
   DHT11_delay_us(2);
   Data=0;
   DHT11_delay_ms(20);   //延时18ms以上
   Data=1;
   DHT11_delay_us(30);
}

uchar DHT11_rec_byte()      //接收一个字节
{
   uchar i,dat=0;
  for(i=0;i<8;i++)    //从高到低依次接收8位数据
   {          
      while(!Data);   ////等待50us低电平过去
      DHT11_delay_us(8);     //延时60us，如果还为高则数据为1，否则为0 
      dat<<=1;           //移位使正确接收8位数据，数据为0时直接移位
      if(Data==1)    //数据为1时，使dat加1来接收数据1
         dat+=1;
      while(Data);  //等待数据线拉低    
    }  
    return dat;
}

void DHT11_receive()      //接收40位的数据
{
    uchar R_H,R_L,T_H,T_L,RH,RL,TH,TL,revise; 
    DHT11_start(); // 向DHT11发送开始信号
    if(Data==0) // DHT11先拉低后拉高总线作为响应
    {
        while(Data==0);   //等待拉高     
        DHT11_delay_us(40);  //拉高后延时80us
        R_H=DHT11_rec_byte();    //接收湿度高八位  
        R_L=DHT11_rec_byte();    //接收湿度低八位  
        T_H=DHT11_rec_byte();    //接收温度高八位  
        T_L=DHT11_rec_byte();    //接收温度低八位
        revise=DHT11_rec_byte(); //接收校正位

        DHT11_delay_us(25);    //结束

        if((R_H+R_L+T_H+T_L)==revise)      //校正
        {
            RH=R_H;
            RL=R_L;
            TH=T_H;
            TL=T_L;
        } 
        /*数据处理，方便显示*/
        rec_dat[0]='l';
        rec_dat[1]='=';
        rec_dat[2]='0'+(RH/10);
        rec_dat[3]='0'+(RH%10);
        rec_dat[4]=',';
        rec_dat[5]='x';
        rec_dat[6]='=';
        rec_dat[7]='0'+(TH/10);
        rec_dat[8]='0'+(TH%10);
        rec_dat[9]='\r';
        rec_dat[10] = '\n';
        rec_dat[11] = '0';
    }
}

void main()
{
   UART_init(); // 串口初始化
   DHT11_delay_ms(1500);    //DHT11上电后要等待1S以越过不稳定状态在此期间不能发送任何指令
   while(1)
   {
       DHT11_receive();
       UART_send_string(rec_dat); //串口发送数据
       delay(); // 延时一定时间
   }
}
```

### UART传输设置

```c
void UART_init()
{
    // 波特率9600
    SCON = 0X50; //工作于方式1  8位无校验异步通信的收发模式，并清空收发中断标志位
    TMOD = 0X20; // 定时器1 8位自动加载计数器
    // 定时器初值，与波特率有关
    TH1 = 0XFD; 
    TL1 = 0XFD; 

    TR1 = 1; // 启动定时器1
}

void UART_send_byte(uchar dat) // 发送一个字节
{
    SBUF = dat;
    while (TI == 0); // 等待直到发送成功
    TI = 0;
}

void UART_send_string(uchar *buf) // 发送字符串
{
    while (*buf != '0')
    {
        UART_send_byte(*buf++);
    }
}
```


