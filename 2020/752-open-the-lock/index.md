# #752 Open the Lock


#### 题目

你有一个带有四个圆形拨轮的转盘锁。每个拨轮都有10个数字： '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' 。每个拨轮可以自由旋转：例如把 '9' 变为  '0'，'0' 变为 '9' 。每次旋转都只能旋转一个拨轮的一位数字。

锁的初始数字为 '0000' ，一个代表四个拨轮的数字的字符串。

列表 deadends 包含了一组死亡数字，一旦拨轮的数字和列表里的任何一个元素相同，这个锁将会被永久锁定，无法再被旋转。

字符串 target 代表可以解锁的数字，你需要给出最小的旋转次数，如果无论如何不能解锁，返回 -1。

```
输入：deadends = ["0201","0101","0102","1212","2002"], target = "0202"
输出：6
解释：
可能的移动序列为 "0000" -> "1000" -> "1100" -> "1200" -> "1201" -> "1202" -> "0202"。
注意 "0000" -> "0001" -> "0002" -> "0102" -> "0202" 这样的序列是不能解锁的，
因为当拨动到 "0102" 时这个锁就会被锁定。
```

提示：

1. 死亡列表 deadends 的长度范围为 [1, 500]。
2. 目标数字 target 不会在 deadends 之中。
3. 每个 deadends 和 target 中的字符串的数字会在 10,000 个可能的情况 '0000' 到 '9999' 中产生。

---

#### 解法1.BFS

直接BFS，不过要加一个集合来表示deadends和重复的pattern。使用`unordered_set<string>`来存储和查询。

```c++
class Solution {
public:
    string add(string s, int i)
    {
        if (s[i] == '9')
            s[i] = '0';
        else
            s[i] = s[i]+1;
        return s;
    }
    
    string minus(string s, int i)
    {
        if (s[i] == '0')
            s[i] = '9';
        else
            s[i] = s[i]-1;
        return s;
    }
    
    int openLock(vector<string>& deadends, string target) {
        unordered_set<string> set(deadends.begin(),deadends.end());
        if (set.find("0000")!=set.end())
            return -1;
        queue<string> q;
        q.push("0000");
        int step = 0;
        
        while (!q.empty())
        {
            int n = q.size();
            for (int j=0; j<n; j++)
            {
                string front = q.front();
                for (int i=0; i<4; i++)
                {
                    string a = add(front,i);
                    if (set.find(a)!=set.end())
                        continue;
                    if (a == target)
                        return step+1;
                    if (!set.count(a))
                    {
                        q.push(a);
                        set.insert(a);
                    }
                }
                for (int i=0; i<4; i++) {
                    string b = minus(front,i);
                    if (set.find(b)!=set.end())
                        continue;
                    if (b == target)
                        return step+1;
                    if (!set.count(b))
                    {
                        q.push(b);
                        set.insert(b);
                    }
                }
                q.pop();
            }
            step+=1;
        }
        return -1;
    }
};
```
#### 解法2.双向BFS


