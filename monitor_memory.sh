#!/bin/bash

# 使用方法函数
usage() {
    echo "使用方法: $0 <程序路径> <输出文件路径> [参数...]"
    echo "例如: $0 ./myprogram /tmp/memory_log.csv arg1 arg2"
    exit 1
}

# 检查参数
if [ $# -lt 2 ]; then
    usage
fi

PROGRAM_PATH="$1"
OUTPUT_FILE="$2"
shift 2  # 移除前两个参数，剩余的是程序参数

# 检查程序是否存在
if [ ! -f "$PROGRAM_PATH" ] || [ ! -x "$PROGRAM_PATH" ]; then
    echo "错误: 程序路径 '$PROGRAM_PATH' 不存在或没有执行权限"
    exit 1
fi

# 检查输出文件路径
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "创建输出目录: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR" || { echo "无法创建输出目录"; exit 1; }
fi

# 初始化输出文件，写入CSV头
echo "时间戳,PID,RSS(KB),VSZ(KB),内存使用率(%)" > "$OUTPUT_FILE"

# 运行程序
echo "启动程序: $PROGRAM_PATH $*"
"$PROGRAM_PATH" "$@" &
PROGRAM_PID=$!

echo "监控PID: $PROGRAM_PID"
echo "内存使用数据将保存到: $OUTPUT_FILE"

# 等待程序启动


# 监控内存使用
SAMPLES=0
while kill -0 $PROGRAM_PID 2>/dev/null; do
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    
    # 使用ps获取RSS(实际物理内存)和VSZ(虚拟内存)
    MEMORY_INFO=$(ps -p $PROGRAM_PID -o pid=,rss=,vsz=,pmem= 2>/dev/null)
    
    if [ -n "$MEMORY_INFO" ]; then
        # 确保正确格式化数据，去除前导空格
        FORMATTED_INFO=$(echo $MEMORY_INFO | tr -s ' ' ',')
        echo "$TIMESTAMP,$FORMATTED_INFO" >> "$OUTPUT_FILE"
        SAMPLES=$((SAMPLES+1))
    else
        # 如果进程已经结束但未被检测到
        break
    fi
    
    sleep 0.001
done

echo "程序已结束，内存使用记录已保存到: $OUTPUT_FILE"
echo "共捕获了 $SAMPLES 个采样点"

# 生成统计信息，但仅在有足够数据时
if [ -f "$OUTPUT_FILE" ] && [ $SAMPLES -gt 0 ]; then
    echo "内存使用统计:"
    
    # 最大RSS
    MAX_RSS=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f3 | sort -n | tail -1)
    echo "最大RSS(KB): $MAX_RSS"
    
    # 平均RSS - 使用awk确保正确计算
    AVG_RSS=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f3 | awk '{ sum += $1; n++ } END { if (n > 0) printf "%.2f", sum / n; else print "0"; }')
    echo "平均RSS(KB): $AVG_RSS"
    
    # 最大内存使用率 - 确保是百分比形式
    MAX_MEM_PERCENT=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f5 | sort -n | tail -1)
    echo "最大内存使用率(%): $MAX_MEM_PERCENT"
    
    # 添加额外信息：最小值和标准差
    MIN_RSS=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f3 | sort -n | head -1)
    echo "最小RSS(KB): $MIN_RSS"
    
    # 内存增长率
    echo "内存增长(KB): $((MAX_RSS - MIN_RSS))"
fi

exit 0