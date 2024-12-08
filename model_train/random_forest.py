import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, auc
import joblib
import os

# 1. 加载数据
primers = pd.read_csv('data/primers_scoring_output.csv')
nonsense_primers = pd.read_csv('data/nonsense_primers_10w_scoring_output.csv')

# 2. 数据预处理：合并数据并标记类别
primers['label'] = 1  # 正常引物标记为1
nonsense_primers['label'] = 0  # 无效引物标记为0

# 3. 合并数据集
data = pd.concat([primers, nonsense_primers])

# 4. 提取特征（从第二列到最后一列）和标签（第一列是lamp_id，最后一列是label）
X = data.iloc[:, 1:-1]  # 特征
y = data['label']  # 标签

# 5. 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 6. 创建并训练随机森林模型
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)

# 7. 在测试集上进行预测
y_pred = rf_model.predict(X_test)
y_pred_proba = rf_model.predict_proba(X_test)[:, 1]  # 获取预测概率

# 8. 输出模型评估结果
print("Classification Report:\n", classification_report(y_test, y_pred))

# 9. 保存训练好的模型
os.makedirs('model', exist_ok=True)
joblib.dump(rf_model, 'model/random_forest_model.pkl')

# 10. 输出特征的重要性（权重）
feature_importances = pd.DataFrame(rf_model.feature_importances_, index=X.columns, columns=["importance"])
print("\nFeature Importances:\n", feature_importances)

# 11. 保存特征重要性到CSV文件
feature_importances.to_csv('model/feature_importances.csv')

# 12. 可视化评估指标
# 混淆矩阵
conf_matrix = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', xticklabels=['Nonsense', 'Primers'], yticklabels=['Nonsense', 'Primers'])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix')
plt.savefig('model/confusion_matrix.pdf')
plt.close()

# ROC曲线
fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
roc_auc = auc(fpr, tpr)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='blue', label=f'ROC Curve (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc='lower right')
plt.savefig('model/roc_curve.pdf')
plt.close()

# 特征重要性
feature_importances = feature_importances.sort_values(by="importance", ascending=False)
plt.figure(figsize=(12, 6))
sns.barplot(x=feature_importances.importance, y=feature_importances.index, palette='viridis', hue=None, legend=False)
plt.xlabel('Feature Importance Score')
plt.ylabel('Features')
plt.title('Feature Importances')
plt.savefig('model/feature_importances.pdf')
plt.close()

print("所有可视化图表已生成并保存到 'model' 文件夹中。")
