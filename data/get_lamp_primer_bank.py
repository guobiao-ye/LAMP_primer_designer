import requests
import pandas as pd
import certifi

# 获取 JSON 数据
url = "https://lampprimer.mathematik.uni-marburg.de/api.php?table=primer&action=list"
response = requests.get(url, verify=False)

data = response.json().get('data', [])  # 使用 get 方法确保即使没有 'data' 也不会报错

# 整理数据
formatted_data = []
for item in data:

    row = {
        "id": item.get("id"),
        "genbank": item.get("genbank"),
        "name": item.get("name"),
        "sequence": item.get("sequence"),
        "position": item.get("position"),
        "lamp_id": item.get("lamp_id")
    }
    formatted_data.append(row)

# 转换为 DataFrame
df = pd.DataFrame(formatted_data)

# 导出为 CSV 文件
df.to_csv("output.csv", index=False, encoding='utf-8-sig')

print("CSV 文件已生成：output.csv")