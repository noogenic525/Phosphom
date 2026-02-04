import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox

def calculate_metrics():
    # 1. GUI 초기화 및 파일 선택
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    
    print("파일을 선택하세요")
    file_path = filedialog.askopenfilename(
        title="선택",
        filetypes=[("Excel files", "*.xlsx *.xls")]
    )
    if not file_path: return

    try:
        # 2. 데이터 로드 (Sheet 2 - 검증 가능 데이터 50개)
        df = pd.read_excel(file_path, sheet_name=1)
        
        # 3. TP, FP, FN 카운트 초기화
        tp = 0 # True Positive
        fp = 0 # False Positive
        fn = 0 # False Negative

        print("\n--- 행별 상세 검증 시작 ---")
        for idx, row in df.iterrows():
            # 예측값 리스트 정제
            p_raw = str(row['Kinase Name']) if pd.notna(row['Kinase Name']) else ""
            p_preds = set([k.strip().upper() for k in p_raw.split(',') if k.strip()])
            
            # 실제 정답 리스트 정제 (Actual_Kinases_in_PSP 컬럼 기준)
            a_raw = str(row['Actual_Kinases_in_PSP']) if pd.notna(row['Actual_Kinases_in_PSP']) else ""
            a_actuals = set([k.strip().upper() for k in a_raw.split(',') if k.strip()])
            
            # 부분 일치(mTOR vs mTOR&CDK)를 고려한 카운팅
            current_row_tps = set()
            
            # TP 찾기: 예측한 것 중 정답에 포함되는 것
            for p in p_preds:
                is_hit = False
                for a in a_actuals:
                    if p in a or a in p:
                        is_hit = True
                        current_row_tps.add(p)
                        break
                if is_hit:
                    tp += 1
                else:
                    fp += 1 # 예측했는데 실제엔 없음
            
            # FN 찾기: 정답 중에 내가 예측하지 못한 것
            for a in a_actuals:
                was_found = False
                for p in p_preds:
                    if p in a or a in p:
                        was_found = True
                        break
                if not was_found:
                    fn += 1 # 정답에 있는데 못 찾음

        # 4. 수식 계산
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

        # 5. 결과 요약
        summary = (
            f"검증 대상: {len(df)}개 모티프\n"
            f"----------------------------\n"
            f"Precision (정밀도): {precision:.4f}\n"
            f"Recall (재현율): {recall:.4f}\n"
            f"F1-score: {f1:.4f}\n"
            f"----------------------------\n"
            f"TP: {tp}, FP: {fp}, FN: {fn}"
        )
        print("\n" + "="*40 + "\n" + summary + "\n" + "="*40)
        messagebox.showinfo("F1-score 검증 완료", summary)

    except Exception as e:
        messagebox.showerror("오류", f"에러가 발생했습니다: {e}")

if __name__ == "__main__":
    calculate_metrics()