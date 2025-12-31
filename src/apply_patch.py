# TK.Labo
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
パッチ適用スクリプト
Usage: python apply_patch.py <patch_file> <target_file> 
"""

import sys
import re
from pathlib import Path

def normalize_line(line):
    """行を正規化（改行コード統一、末尾空白除去）"""
    # 改行コードを統一
    return line.replace('\r\n', '\n').replace('\r', '\n')

def parse_unified_diff(patch_content):
    """Parse unified diff format patch"""
    lines = normalize_line(patch_content).strip().split('\n')
    hunks = []
    current_hunk = None
    
    for line in lines:
        if line.startswith('@@'):
            match = re.match(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@', line)
            if match:
                old_start = int(match.group(1))
                old_count = int(match.group(2)) if match.group(2) else 1
                new_start = int(match.group(3))
                new_count = int(match.group(4)) if match.group(4) else 1
                
                if current_hunk:
                    hunks.append(current_hunk)
                
                current_hunk = {
                    'old_start': old_start,
                    'old_count': old_count,
                    'new_start': new_start,
                    'new_count': new_count,
                    'lines': []
                }
        elif current_hunk is not None:
            if line.startswith('---') or line.startswith('+++'):
                continue
            if line.startswith(' ') or line.startswith('+') or line.startswith('-'):
                current_hunk['lines'].append(line)
            elif line.startswith('\\'):
                continue
    
    if current_hunk:
        hunks.append(current_hunk)
    
    return hunks

def verify_hunk_context(target_lines, hunk, start_pos):
    """ハンクがこの位置に適用可能か検証"""
    target_idx = start_pos
    
    for line in hunk['lines']:
        if line.startswith(' ') or line.startswith('-'):
            # コンテキスト行または削除行は元ファイルと一致すべき
            if target_idx >= len(target_lines):
                return False
            
            expected = normalize_line(line[1:])
            actual = normalize_line(target_lines[target_idx].rstrip('\n'))
            
            # 厳密比較と緩い比較の両方を試す
            if expected.rstrip() != actual.rstrip():
                print(f"  ミスマッチ検出 行{target_idx + 1}:")
                print(f"    期待: '{expected.rstrip()}'")
                print(f"    実際: '{actual.rstrip()}'")
                return False
            
            target_idx += 1
        elif line.startswith('+'):
            # 追加行は検証不要
            continue
    
    return True

def find_hunk_position(target_lines, hunk, suggested_pos):
    """ハンクを適用できる位置を探す"""
    # まず指定位置を試す
    if verify_hunk_context(target_lines, hunk, suggested_pos):
        return suggested_pos
    
    print(f"  警告: 行{suggested_pos + 1}では適用できません。近傍を検索中...")
    
    # 前後10行の範囲で検索
    for offset in range(-10, 11):
        test_pos = suggested_pos + offset
        if 0 <= test_pos < len(target_lines):
            if verify_hunk_context(target_lines, hunk, test_pos):
                print(f"  → 行{test_pos + 1}で適用可能な位置を発見")
                return test_pos
    
    return None

def apply_patch(target_file, patch_file, output_file=None, fuzzy=True):
    """Apply patch to target file"""
    
    # Read target file
    try:
        with open(target_file, 'r', encoding='utf-8') as f:
            target_lines = [normalize_line(line) for line in f.readlines()]
    except UnicodeDecodeError:
        try:
            with open(target_file, 'r', encoding='shift_jis') as f:
                target_lines = [normalize_line(line) for line in f.readlines()]
        except UnicodeDecodeError:
            with open(target_file, 'r', encoding='cp932') as f:
                target_lines = [normalize_line(line) for line in f.readlines()]
    
    # Read patch file
    with open(patch_file, 'r', encoding='utf-8') as f:
        patch_content = f.read()
    
    # Parse patch
    hunks = parse_unified_diff(patch_content)
    
    if not hunks:
        print("エラー: パッチファイルにハンクが見つかりません")
        return False
    
    print(f"パッチ解析完了: {len(hunks)}個のハンクを検出")
    print(f"対象ファイル: {len(target_lines)}行\n")
    
    # Apply hunks
    line_offset = 0
    
    for i, hunk in enumerate(hunks, 1):
        print(f"ハンク {i}/{len(hunks)} 処理中...")
        
        # Calculate adjusted position
        suggested_pos = hunk['old_start'] - 1 + line_offset
        
        # Find correct position (with fuzzy matching if enabled)
        if fuzzy:
            actual_pos = find_hunk_position(target_lines, hunk, suggested_pos)
            if actual_pos is None:
                print(f"エラー: ハンク{i}を適用できる位置が見つかりません")
                print(f"  期待位置: 行{hunk['old_start']}")
                return False
        else:
            if not verify_hunk_context(target_lines, hunk, suggested_pos):
                print(f"エラー: ハンク{i}のコンテキストが一致しません")
                return False
            actual_pos = suggested_pos
        
        # Build new content
        new_content = []
        old_line_idx = actual_pos
        
        for line in hunk['lines']:
            if line.startswith(' '):
                # Context line
                if old_line_idx < len(target_lines):
                    new_content.append(target_lines[old_line_idx])
                    old_line_idx += 1
            elif line.startswith('-'):
                # Delete line
                old_line_idx += 1
            elif line.startswith('+'):
                # Add line
                content = line[1:]
                if not content.endswith('\n'):
                    content += '\n'
                new_content.append(content)
        
        # Calculate replacement range
        lines_to_remove = sum(1 for l in hunk['lines'] if l.startswith(' ') or l.startswith('-'))
        end_idx = actual_pos + lines_to_remove
        
        # Apply hunk
        target_lines[actual_pos:end_idx] = new_content
        
        # Update offset
        lines_added = len(new_content)
        line_offset += (lines_added - lines_to_remove)
        
        print(f"  - 適用完了: 行{actual_pos + 1} ({lines_to_remove}行削除, {lines_added}行追加)\n")
    
    # Determine output file
    if output_file is None:
        # Default: overwrite original (with backup)
        backup_file = target_file + '.orig'
        if Path(backup_file).exists():
            import datetime
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_file = f"{target_file}.orig.{timestamp}"
        
        Path(target_file).rename(backup_file)
        print(f" - バックアップ作成: {backup_file}")
        output_path = target_file
    else:
        # Write to new file
        output_path = output_file
        print(f" - 元のファイルは保持: {target_file}")
    
    # Write result
    with open(output_path, 'w', encoding='utf-8', newline='\n') as f:
        f.writelines(target_lines)
    
    print(f"- パッチ適用完了: {output_path}\n")
    return True

def main():
    if len(sys.argv) < 3:
        print("使い方: python apply_patch.py <patch_file> <target_file> [output_file] [--strict]")
        print("  output_file: 出力ファイル名（省略時は元ファイルを上書き）")
        print("  --strict: あいまい検索を無効化（完全一致のみ）")
        print("\n例:")
        print("  python apply_patch.py modify_rri.patch 1.4.2.7/RRI.f90")
        print("  python apply_patch.py modify_rri.patch 1.4.2.7/RRI.f90 1.4.2.7/RRI_UNST.f90")
        print("  python apply_patch.py modify_rri.patch 1.4.2.7/RRI.f90 1.4.2.7/RRI_UNST.f90 --strict")
        sys.exit(1)
    
    patch_file = sys.argv[1]
    target_file = sys.argv[2]
    
    # Parse optional arguments
    output_file = None
    fuzzy = True
    
    for arg in sys.argv[3:]:
        if arg == '--strict':
            fuzzy = False
        elif not arg.startswith('--'):
            output_file = arg
    
    if not Path(patch_file).exists():
        print(f"エラー: パッチファイル '{patch_file}' が見つかりません")
        sys.exit(1)
    
    if not Path(target_file).exists():
        print(f"エラー: 対象ファイル '{target_file}' が見つかりません")
        sys.exit(1)
    
    try:
        success = apply_patch(target_file, patch_file, output_file, fuzzy)
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\nエラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
