# #739 每日温度


单调栈解法

<!--more-->

#### 题目

请根据每日 `气温` 列表，重新生成一个列表。对应位置的输出为：要想观测到更高的气温，至少需要等待的天数。如果气温在这之后都不会升高，请在该位置用 `0` 来代替。

例如，给定一个列表 `temperatures = [73, 74, 75, 71, 69, 72, 76, 73]`，你的输出应该是 `[1, 1, 4, 2, 1, 1, 0, 0]`。

提示：气温 列表长度的范围是 `[1, 30000]`。每个气温的值的均为华氏度，都是在 `[30, 100]` 范围内的整数。

#### 单调栈解法

维护一个栈，使得从栈底到栈顶始终是递增的，具体做法举例：

假如列表`temperatures = [10,12,5,7,13,14]`

- 首先10入栈，12比10大，10弹出，栈顶为12 `[12]`

- 5比12小，入栈 `[12 5]`
- 7比5大，5弹出，7入栈 `[12 7]`
- 13比7大，7弹出，13比12大，12弹出，13入栈 `[13]`
- 14比13大，13弹出，14入栈  `[14]`

为了获取等待天数，在每个元素出栈时计算此元素与与之比较的元素的下标差。为简单起见，栈中存储下标而不是值。

所以实际过程如下：

{{< figure src="/images/739-daily-temperatures/stack.jpg" title="操作顺序" >}}

代码：

```c++
class Solution {
public:
    vector<int> dailyTemperatures(vector<int>& T) {
        stack<int> s;
        int n = T.size();
        // initionlize
        vector<int> ans(n,0);
        
        for (int i=0; i< n; ++i) {
            while(!s.empty() && T[i] > T[s.top()]) {
                int index = s.top();
                // inside s are indexes
                ans[index] = (i-s.top());
                s.pop();
            }
            // got a small number, then push it to stack
            s.push(i);
        }
        return ans;
    }
};
```

#### 复杂度分析

- 时间复杂度：$O(N)$ ，温度列表的所有元素都将进栈一次
- 空间复杂度:  $O(N)$，结果列表的长度为$N$
