#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple patch file applicator for Windows
Usage: python apply_patch.py <patch_file> <target_file>
Example: python apply_patch.py modify_rri.patch RRI.f90
"""

import sys
import re
from pathlib import Path

def parse_unified_diff(patch_content):
    """Parse unified diff format patch"""
    lines = patch_content.strip().split('\n')
    hunks = []
    current_hunk = None
    
    for line in lines:
        # Hunk header: @@ -start,count +start,count @@
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
            # Skip file headers
            if line.startswith('---') or line.startswith('+++'):
                continue
            # Context, addition, or deletion
            if line.startswith(' ') or line.startswith('+') or line.startswith('-'):
                current_hunk['lines'].append(line)
            elif line.startswith('\\'):
                # "\ No newline at end of file"
                continue
    
    if current_hunk:
        hunks.append(current_hunk)
    
    return hunks

def apply_patch(target_file, patch_file):
    """Apply patch to target file"""
    
    # Read files
    try:
        with open(target_file, 'r', encoding='utf-8') as f:
            target_lines = f.readlines()
    except UnicodeDecodeError:
        # Try with different encoding
        with open(target_file, 'r', encoding='shift_jis') as f:
            target_lines = f.readlines()
    
    with open(patch_file, 'r', encoding='utf-8') as f:
        patch_content = f.read()
    
    # Parse patch
    hunks = parse_unified_diff(patch_content)
    
    if not hunks:
        print("警告: パッチファイルにハンク(変更箇所)が見つかりませんでした")
        print("パッチファイルの形式を確認してください")
        return False
    
    # Apply hunks with cumulative offset
    line_offset = 0
    
    for hunk in hunks:
        # Adjust starting position based on previous changes
        old_start = hunk['old_start'] - 1 + line_offset  # Convert to 0-based
        
        # Build new content for this hunk
        new_content = []
        old_line_idx = old_start
        
        for line in hunk['lines']:
            if line.startswith(' '):
                # Context line - copy from original
                if old_line_idx < len(target_lines):
                    new_content.append(target_lines[old_line_idx])
                    old_line_idx += 1
            elif line.startswith('-'):
                # Line to remove - skip it
                old_line_idx += 1
            elif line.startswith('+'):
                # Line to add
                content = line[1:]
                if not content.endswith('\n'):
                    content += '\n'
                new_content.append(content)
        
        # Calculate how many lines to replace
        lines_to_remove = sum(1 for l in hunk['lines'] if l.startswith(' ') or l.startswith('-'))
        
        # Replace the section
        end_idx = old_start + lines_to_remove
        target_lines[old_start:end_idx] = new_content
        
        # Update offset for next hunk
        lines_added = len(new_content)
        line_offset += (lines_added - lines_to_remove)
        
        print(f"ハンク適用: 行 {hunk['old_start']} 付近 ({lines_to_remove}行削除, {lines_added}行追加)")
    
    # Write back
    backup_file = target_file + '.orig'
    if Path(backup_file).exists():
        # If backup exists, add timestamp
        import datetime
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_file = f"{target_file}.orig.{timestamp}"
    
    Path(target_file).rename(backup_file)
    print(f"元のファイルを {backup_file} にバックアップしました")
    
    with open(target_file, 'w', encoding='utf-8') as f:
        f.writelines(target_lines)
    
    print(f"パッチを {target_file} に適用しました")
    return True

def main():
    if len(sys.argv) != 3:
        print("使い方: python apply_patch.py <patch_file> <target_file>")
        print("例: python apply_patch.py modify_rri.patch RRI.f90")
        sys.exit(1)
    
    patch_file = sys.argv[1]
    target_file = sys.argv[2]
    
    if not Path(patch_file).exists():
        print(f"エラー: パッチファイル '{patch_file}' が見つかりません")
        sys.exit(1)
    
    if not Path(target_file).exists():
        print(f"エラー: 対象ファイル '{target_file}' が見つかりません")
        sys.exit(1)
    
    try:
        apply_patch(target_file, patch_file)
    except Exception as e:
        print(f"エラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()