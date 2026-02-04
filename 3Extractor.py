import re
import pandas as pd
from docx import Document
from tkinter import filedialog, messagebox
import tkinter as tk

# 사용자 정의 Kinase 키워드 (정답지 시트 및 분석 로직과 일치)
KINASE_KEYS = [
    "IKK", "GSK-3beta&CDK", "GSK-3beta", "Erk", "AMPK", 
    "p38&MAPK", "ATM&ATR", "PKA", "CaMK2", "CK2", "mTOR&CDK"
]

def get_clean_gene_symbol(text):
    """긴 헤더 문장에서 (Syn1) 같은 유전자 심볼만 추출"""
    # 괄호 안의 내용을 찾되, 숫자로만 이루어지지 않은 것(심볼)을 우선 탐색
    match = re.search(r'\(([^0-9][^)]+)\)', text)
    if match:
        return match.group(1).strip()
    # 괄호가 없다면 첫 단어 혹은 전체 텍스트 반환
    return text.split(',')[0].strip()

def clean_motif_sequence(text):
    """서열 내의 (0.477) 같은 확률값과 공백을 제거"""
    # 괄호 안의 소숫점 숫자 패턴 제거
    text = re.sub(r'\(0\.\d+\)', '', text)
    text = re.sub(r'\(1\)', '', text)
    # 아미노산 대문자 제외한 나머지(숫자, 기호, 공백) 제거
    text = re.sub(r'[^A-Z]', '', text.upper())
    return text.strip()

def extract_kinases_from_info(info_text):
    """문장 끝 괄호 (765 CK2, 766 CK2) 에서 Kinase 이름만 추출"""
    found = []
    upper_info = info_text.upper()
    for key in KINASE_KEYS:
        # 키워드가 텍스트 안에 포함되어 있는지 확인
        if key.upper() in upper_info:
            found.append(key)
    # 중복 제거 및 정렬
    return ", ".join(dict.fromkeys(found))

def is_header_line(text):
    """유전자 헤더 라인인지 판별 (보통 mRNA, variant 등을 포함함)"""
    keywords = ['mRNA', 'variant', 'complete cds', 'protein', 'Syn1', 'Syngap1', 'Septin7']
    return any(k.lower() in text.lower() for k in keywords)

def run_extractor():
    root = tk.Tk()
    root.withdraw()
    
    word_file = filedialog.askopenfilename(title="분석된 Word 파일 선택", filetypes=[("Word", "*.docx")])
    if not word_file: return
    
    save_file = filedialog.asksaveasfilename(title="엑셀 저장 경로 선택", defaultextension=".xlsx", filetypes=[("Excel", "*.xlsx")])
    if not save_file: return

    try:
        doc = Document(word_file)
        extracted_data = []
        current_gene_symbol = "Unknown"

        print("데이터 정밀 분석 중...")

        for para in doc.paragraphs:
            text = para.text.strip()
            if not text: continue

            # 1. 유전자 헤더 라인인 경우 심볼 업데이트
            if is_header_line(text) and len(text) > 20:
                current_gene_symbol = get_clean_gene_symbol(text)
                continue

            # 2. 모티프 라인 판별
            # 특징: 문장 끝이 )로 끝나고, 내부에 숫자+Kinase 이름이 있거나 (0.xxx)가 있음
            if "(" in text and ")" in text:
                # 괄호가 여러 개일 수 있으므로 마지막 괄호를 정보창으로 분리
                parts = text.rsplit("(", 1)
                motif_part = parts[0].strip()
                info_part = parts[1].strip()

                # 서열 정제
                clean_seq = clean_motif_sequence(motif_part)
                
                # Kinase 추출
                kinases = extract_kinases_from_info(info_part)
                
                # 전체 서열 블록(매우 긴 것)은 제외 (보통 모티프는 50자 미만)
                if clean_seq and 3 < len(clean_seq) < 100:
                    extracted_data.append({
                        "Gene Name": current_gene_symbol,
                        "Motif": clean_seq,
                        "Kinase Name": kinases
                    })

        # 데이터프레임 생성 및 저장
        df = pd.DataFrame(extracted_data)
        
        if df.empty:
            messagebox.showwarning("알림", "조건에 맞는 모티프 데이터를 찾지 못했습니다.")
            return

        df.to_excel(save_file, index=False)
        
        summary = f"추출 완료!\n유전자 수: {df['Gene Name'].nunique()}개\n모티프 행 수: {len(df)}개"
        print(summary)
        messagebox.showinfo("성공", summary)

    except Exception as e:
        messagebox.showerror("오류", f"에러 발생: {e}")

if __name__ == "__main__":
    run_extractor()