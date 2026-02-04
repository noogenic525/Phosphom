import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox
import os

def run_validator_v3():
    # 1. GUI 초기화 및 파일 선택
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    
    print("추출된 내 데이터 파일을 선택하세요...")
    my_file_path = filedialog.askopenfilename(
        title="내 데이터 선택",
        filetypes=[("Excel files", "*.xlsx *.xls")]
    )
    if not my_file_path: return

    print("정답(PSP) 엑셀 파일을 선택하세요...")
    ref_file_path = filedialog.askopenfilename(
        title="정답 엑셀 선택 (Substrates of protein)",
        filetypes=[("Excel files", "*.xlsx *.xls")]
    )
    if not ref_file_path: return

    try:
        # 2. 내 데이터 로드
        pm_df = pd.read_excel(my_file_path)
        
        # 3. 정답지 시트 로드 및 데이터 구조화
        # 시트 이름과 매칭할 키워드들 (내 DB의 이름들)
        kinase_keywords = [
            "IKK", "GSK-3beta&CDK", "GSK-3beta", "Erk", "AMPK", 
            "p38&MAPK", "ATM&ATR", "PKA", "CaMK2", "CK2", "mTOR&CDK"
        ]

        excel_reader = pd.ExcelFile(ref_file_path)
        all_sheets = excel_reader.sheet_names
        
        # kinase_ref_data['PKA'] = [서열1, 서열2, ...] 형태로 저장
        kinase_ref_data = {}
        
        print("\n--- 정답지 시트 스캔 및 로딩 ---")
        for k_name in kinase_keywords:
            # 시트 이름 중 키워드가 포함된 시트 찾기 (대소문자 무시)
            # 예: ATM&ATR 키워드로 AMT&ATR 시트도 찾을 수 있게 유연하게 처리
            target_sheet = None
            for s in all_sheets:
                # & 기호나 철자 오타를 고려하여 키워드의 주요 부분이 포함되는지 확인
                clean_k = k_name.replace('&', '').upper()
                clean_s = s.replace('&', '').upper()
                if clean_k in clean_s or clean_s in clean_k:
                    target_sheet = s
                    break
            
            if target_sheet:
                df_ref = pd.read_excel(ref_file_path, sheet_name=target_sheet)
                if 'SITE_+/-7_AA' in df_ref.columns:
                    # 언더바 제거 및 대문자화하여 서열 리스트 생성
                    seqs = df_ref['SITE_+/-7_AA'].astype(str).str.replace('_', '').str.upper().unique().tolist()
                    kinase_ref_data[k_name] = seqs
                    print(f"✅ [{k_name}] 매칭 성공 -> 시트명: {target_sheet} ({len(seqs)}개 서열)")
                else:
                    print(f"❌ [{target_sheet}] 시트에 'SITE_+/-7_AA' 컬럼이 없습니다.")
        
        # 4. 검증 로직 수행
        results = []
        print("\n--- 데이터 검증 시작 ---")
        
        for idx, row in pm_df.iterrows():
            # 내 데이터의 모티프와 예측값
            p_motif = str(row['Motif']).upper()
            p_preds_raw = str(row['Kinase Name']) if pd.notna(row['Kinase Name']) else ""
            p_preds = [p.strip().upper() for p in p_preds_raw.split(',') if p.strip()]
            
            # 이 모티프가 어떤 정답 시트들에 들어있는지 전수 조사
            found_kinases = []
            for k_name, ref_seqs in kinase_ref_data.items():
                # 포함 관계 확인 (내 모티프가 정답 서열에 있거나, 정답 서열이 내 모티프에 있거나)
                if any(p_motif in s or s in p_motif for s in ref_seqs):
                    found_kinases.append(k_name)
            
            # 정답 판정 (부분 일치 로직 적용)
            is_correct = False
            if found_kinases:
                actual_upper_list = [a.upper() for a in found_kinases]
                for pred in p_preds:
                    if is_correct: break
                    for actual in actual_upper_list:
                        if pred in actual or actual in pred:
                            is_correct = True
                            break
            
            results.append({
                'In_PSP': len(found_kinases) > 0,
                'Correct': is_correct,
                'Actual_Kinases_in_PSP': ", ".join(found_kinases)
            })

        # 5. 결과 합치기 및 정확도 계산
        res_df = pd.DataFrame(results)
        final_df = pd.concat([pm_df, res_df], axis=1)
        
        valid_rows = res_df[res_df['In_PSP'] == True]
        acc = (valid_rows['Correct'].sum() / len(valid_rows) * 100) if not valid_rows.empty else 0
        
        # 6. 결과 저장
        save_path = filedialog.asksaveasfilename(
            title="검증 결과 저장",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx")]
        )
        
        if save_path:
            final_df.to_excel(save_path, index=False)
            
            summary = (
                f"전체 모티프: {len(pm_df)}개\n"
                f"정답지(PSP)와 매칭된 모티프: {len(valid_rows)}개\n"
                f"예측 성공(Correct): {valid_rows['Correct'].sum()}개\n"
                f"최종 정확도(ACC): {acc:.2f}%"
            )
            print("\n" + "="*40 + "\n" + summary + "\n" + "="*40)
            messagebox.showinfo("검증 완료", summary)

    except Exception as e:
        import traceback
        print(traceback.format_exc())
        messagebox.showerror("오류", f"실행 중 문제가 발생했습니다:\n{e}")

if __name__ == "__main__":
    run_validator_v3()