import csv

# 输入文件名和输出文件名
input_file = "output.csv"
output_file = "restructured_output.csv"

def parse_position_and_sequence(name, position, sequence):
    """
    根据位置拆解序列
    :param position: 形如"2313-2337+TTTT+2263-2284"的字符串
    :param sequence: 全部序列
    :return: [(标签1, 序列片段1, 位置1), (标签2, 序列片段2, 位置2)]
    """
    # 处理位置部分，跳过"TTTT"
    parts = position.split("+")
    segments = []

    current_index = 0
    for i, part in enumerate(parts):
        if part == "TTTT":
            continue  # 跳过"TTTT"
        if part.strip():  # 检查是否为空
            try:
                start, end = map(int, part.split("-"))
                length = end - start + 1
                segments.append((f"F{i+1}" if "FIP" == name else f"B{i+1}", sequence[current_index:current_index+length], part))
                current_index += length
            except ValueError:
                print(f"Warning: Invalid position format: {part}")
                continue  # 如果格式错误，跳过该部分
    return segments

def process_csv(input_file, output_file):
    with open(input_file, "r", encoding="utf-8-sig") as infile, open(output_file, "w", newline="", encoding="utf-8") as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ["id", "genbank", "name", "sequence", "position", "lamp_id"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        writer.writeheader()

        for row in reader:
            # Check if any column in the row is empty, and skip the row if so
            if any(value.strip() == "" for value in row.values()):
                continue  # Skip this row if any field is empty

            if row["name"] in ["FIP", "BIP"]:
                segments = parse_position_and_sequence(row['name'], row["position"], row["sequence"])
                for label, seq, pos in segments:
                    new_row = {
                        "id": row["id"],
                        "genbank": row["genbank"],
                        "name": label,
                        "sequence": seq,
                        "position": pos,  # 对应的位置
                        "lamp_id": row["lamp_id"]
                    }
                    writer.writerow(new_row)

            elif row["name"] in ["F3", "B3"]:
                # 保留F3和B3行
                writer.writerow(row)

if __name__ == "__main__":
    process_csv(input_file, output_file)
    print(f"文件已重构并保存为 {output_file}")
