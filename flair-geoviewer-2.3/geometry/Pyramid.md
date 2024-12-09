### 平行于坐标轴的棱锥（Pyramid parallel to a coordinate axis）描述内容的详细解释：

### 基本定义
- **类型说明**： `PYX`、`PYY`、`PYZ` 分别代表平行于 `x` 轴、`y` 轴、`z` 轴的棱锥类型。例如 `PYX` 类型的棱锥，其形状特点是棱锥的轴线方向与 `x` 轴平行，`PYY` 和 `PYZ` 同理，只是各自对应的平行坐标轴不同。
- **形态情况**：这类棱锥可以是完整的（full），也可以是被截断的（truncated）形态，即可能存在顶部被削去一部分的情况，类似棱台的样子。

### 参数含义
- **中心坐标参数**：
    - `V_x`、`V_y`、`V_z` 这三个参数用于确定棱锥大底面（major base）的中心坐标位置。比如对于一个 `PYZ` 类型的棱锥，给出的 `V_x`、`V_y`、`V_z` 值就明确了其大底面中心在三维空间中的具体位置，无论这个棱锥是完整的还是截断的，大底面中心坐标都是由这三个参数来定义的。
- **底面边长参数**：
    - 对于不同类型的棱锥有不同的边长参数组合。以 `PYZ` 为例，`L_y`、`L_z` 分别表示大底面沿着 `y` 方向和 `z` 方向的半边长（half-lengths）。也就是说，如果要确定大底面在 `y` 方向的实际长度范围，需要将 `L_y` 的值乘以 2，`z` 方向同理。对于 `PYX` 类型，则是 `L_z`、`L_x` 来分别表示相应方向的半边长，`PYY` 类型则是 `L_x`、`L_y`。通过这些参数可以准确描述出大底面在对应坐标轴平面内的尺寸大小。
- **高度参数**：
    - `H` 代表棱锥的高度，其数值的正负还能体现方向，具体是从大底面（major basis）指向小底面（minor basis）的方向。正的 `H` 值表示按照规定的正向（比如沿坐标轴正方向）从小底面到大底面，负的则相反，通过这个符号约定能清晰表示棱锥的空间姿态。
- **截断相关参数**：
    - `R` 是一个表示比例的参数，用于描述小底面的边长与大底面边长之间的比例关系。如果 `R = 0.0`，意味着这个棱锥没有被截断，就是一个完整的棱锥；而当 `R >= 1` 时是不符合规定的情况，不被接受，因为从常规理解来说，小底面边长如果大于等于大底面边长就不符合棱锥截断的合理形态了。

### 示例说明
- 例如给出的示例 `PYZ body3        0.0   0.5   -10.0   +20.0  +40.0  50.0 0.8`：
    - 首先表明这是一个 `PYZ` 类型的棱锥，名称为 `body3`。
    - 其大底面中心坐标为 `(0.0, 0.5, -10.0)`，意味着在空间中，大底面中心处于 `x = 0.0`、`y = 0.5`、`z = -10.0` 的位置。
    - 大底面沿着 `x` 方向的半边长 `L_y` 对应的实际范围是从 `x = -20.0` 到 `x = +20.0`（因为 `L_y` 值为 `+20.0`，实际长度是 `2×L_y`），沿着 `y` 方向的半边长 `L_z` 对应的实际范围是从 `y = -39.5` 到 `y = 40.5`（考虑到中心 `y` 坐标是 `0.5`，两边各延伸 `L_z` 值 `40.0`，所以范围是 `0.5 - 40.0` 到 `0.5 + 40.0`）。
    - 棱锥的高度 `H` 为 `50.0`（单位为厘米，文中默认的单位，当然具体单位也可根据实际应用场景判断），说明沿 `z` 轴方向从小底面到大底面的距离是 `50` 厘米。
    - 截断比例 `R` 为 `0.8`，表明这是一个截断的棱锥，小底面的边长与大底面边长存在 `0.8` 的比例关系，进而可以推算出小底面在 `x` 方向延伸范围是从 `-16.0` 到 `+16.0`（是大底面 `x` 方向边长范围 `(-20.0` 到 `+20.0)` 的 `0.8` 倍），在 `y` 方向延伸范围是从 `-31.5` 到 `32.5`（同样基于与大底面 `y` 方向边长的比例关系计算得出）。

### 英文：

Pyramid parallel to a coordinate axis. Code: PYX, PYY, PYZ
-      A PYX (PYY, PYZ)  is a pyramid  parallel to the x (y, z) axis.
-      It Can be full or truncated
-      Each  PYX (PYY, PYZ)  is defined by 7 parameters:
      V_x, V_y, V_z, (coordinates of the centre of the major base)
-      L_y, L_z (res.p L_z, L_x or L_x, L_y) half-lengths of the sides of the major base
-      H : height of the pyramid. The sign gives the direction from major to minor basis
-      R : ratio between the sides of the minor base and the sides of the major base. If R=0.0, the pyramid is not truncated. R>=1.0 is not accepted
-      Example
-      PYZ body3        0.0   0.5   -10.0   +20.0  +40.0  50.0 0.8
-      Is a truncated pyramid with  major base centered at 0.0, 0.5,-10.0 and extending from x=-20.0 to x=+20, and from y=-39.5 to y=40.5. Its axis is along z, its height is 50 cm, the minor basis extends from -16.0 to +16.0 in x and from -31.5 to 32.5 in y