# #200 岛屿数量

最佳解法：广度优先搜素
<!--more-->

> 给你一个由 `'1'`（陆地）和 `'0'`（水）组成的的二维网格，请你计算网格中岛屿的数量。
>
> 岛屿总是被水包围，并且每座岛屿只能由水平方向或竖直方向上相邻的陆地连接形成。
>
> 此外，你可以假设该网格的四条边均被水包围。

示例:

```
输入:
11110
11010
11000
00000
输出: 1
```

```
输入:
11000
11000
00100
00011
输出: 3
解释: 每座岛屿只能由水平和/或竖直方向上相邻的陆地连接而成。
```

**思路**：

- 扫描一遍整个网格，遇到陆地就对他广度优先搜索，并将搜索到的全部陆地置为`0`，这样执行搜索的次数便是岛屿的数量。

- 使用队列实现广搜，按照距离根节点距离递增的顺序遍历节点，遇到陆地则该点坐标入队，同时其上一个节点坐标出队，当队列为空时搜索完成。
- 坐标点表示使用STL容器`pair<int,int>`，头文件`utility`
- 时间复杂度$O(MN)$，扫描二维网格的时间复杂度
- 最坏情况为全是陆地，因此空间复杂度$O(min(M,N))$

**实现**

```c++
class Solution {
public:
    int numIslands(vector<vector<char>>& grid) {
        int nr = grid.size();
        // 经典陷阱
        if (nr == 0)
            return 0;
        int nc = grid[0].size();
        int num = 0;

        for (int r=0; r < nr; r++) {
            for (int c=0; c < nc; c++) {
                if (grid[r][c] == '1') {
                    // 遇到根节点为陆地，岛屿加一
                    num++;
                    grid[r][c] = '0';
                    // 节点坐标
                    queue<pair<int,int>> q;
                    q.push({r,c});
                    // 广搜
                    while(!q.empty()) {
                        auto rc = q.front();
                        int row = rc.first, col = rc.second;
                        // 访问四个邻接点
                        if (row+1 < nr && grid[row+1][col] == '1') {
                            grid[row+1][col] = '0';
                            q.push({row+1,col});
                        }
                        if (row-1 >= 0 && grid[row-1][col] == '1') {
                            grid[row-1][col] = '0';
                            q.push({row-1,col});
                        }
                        if (col+1 < nc && grid[row][col+1] == '1') {
                            grid[row][col+1] = '0';
                            q.push({row,col+1});
                        }
                        if (col-1 >= 0 && grid[row][col-1] == '1') {
                            grid[row][col-1] = '0';
                            q.push({row,col-1});
                        }
                        q.pop();
                    }
                }
            }
        }
        return num;
    }
};
```

**总结**

广搜的实现还是较为容易，第一次做时以为只需要向根节点右下方遍历就行了，结果遇到了特例：*工字形* 的岛屿等，如

```
111
010
111
```

若是只向右下方搜索就会得出岛屿为2。


