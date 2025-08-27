from playwright.sync_api import sync_playwright
import json
import csv
import re

# Output CSV file
output_file = "acp_ks_pks_only.csv"

# Step 1: Load BGC IDs of type 'pks' from local JSON file
def get_pks_bgc_ids_from_file(filepath):
    print("[i] Reading local MIBiG JSON file...")
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    pks_ids = []
    for entry in data["records"]:
        types = entry.get("biosyn_class", [])
        if "pks" in [t.lower() for t in types]:
            pks_ids.append(entry["mibig_accession"])

    print(f"[✓] Found {len(pks_ids)} PKS-only BGCs")
    return pks_ids

# Step 2: Extract tooltip data
def extract_info_from_tooltip(html):
    if re.search(r'^PKS_KS', html):
        domain_type = 'KS'
    elif re.search(r'^ACP', html) or re.search(r'^PKS_PP', html) or re.search(r'^PP-binding', html) or re.search(r'^ACP_beta', html):
        domain_type = 'ACP'
    else:
        return None

    location_match = re.search(r'Location:\s*([\d\-]+ AA)', html)
    aa_match = re.search(r'AA sequence:.*?data-seq="([^"]+)"', html)
    nt_match = re.search(r'Nucleotide sequence:.*?data-seq="([^"]+)"', html)

    return (
        domain_type,
        location_match.group(1) if location_match else "Unknown",
        aa_match.group(1) if aa_match else "Not found",
        nt_match.group(1) if nt_match else "Not found"
    )

# Step 3: Process each BGC
def process_bgc(page, bgc_id):
    url = f"https://mibig.secondarymetabolites.org/repository/{bgc_id}"
    print(f"\n[→] Processing {url}")
    page.goto(url, timeout=60000)
    page.wait_for_timeout(3000)

    try:
        page.get_by_text("NRPS/PKS domains").click(timeout=5000)
        page.wait_for_timeout(2000)
    except:
        print(f"[!] Could not open domain view for {bgc_id}")
        return []

    page.evaluate("""
        document.querySelectorAll('div.jsdomain-tooltip').forEach(tip => {
            tip.style.display = 'block';
            tip.style.visibility = 'visible';
            tip.style.opacity = 1;
        });
    """)

    gene_data = []
    for label in page.locator("text.jsdomain-orflabel").all():
        name = label.text_content().strip()
        box = label.bounding_box()
        if box:
            gene_data.append({
                "name": name,
                "x": box["x"],
                "y": box["y"],
                "y_mid": box["y"] + box["height"] / 2
            })

    module_rects = []
    module_counter = 1
    for mod in page.locator("rect.jsdomain-module, rect.jsdomain-incomplete-module").all():
        box = mod.bounding_box()
        if box:
            box["module_number"] = module_counter
            module_rects.append(box)
            module_counter += 1

    tooltips = page.locator("div.jsdomain-tooltip")
    rects = page.locator("rect.jsdomain-domain")

    domain_entries = []

    for i in range(tooltips.count()):
        html = tooltips.nth(i).inner_html()
        result = extract_info_from_tooltip(html)
        if not result:
            continue
        domain_type, location, aa_seq, nt_seq = result
        box = rects.nth(i).bounding_box()
        if not box:
            continue
        x, y = box["x"], box["y"]
        y_mid = y + box["height"] / 2

        gene_name = "Unknown"
        for g in gene_data:
            if abs(g["y_mid"] - y_mid) < 10:
                gene_name = g["name"]
                break

        module = 0
        for mod in module_rects:
            mod_x1, mod_x2 = mod["x"], mod["x"] + mod["width"]
            mod_y1, mod_y2 = mod["y"], mod["y"] + mod["height"]
            if (x >= mod_x1 and x <= mod_x2) and (y >= mod_y1 and y <= mod_y2):
                module = mod["module_number"]
                break

        domain_entries.append({
            "BGC ID": bgc_id,
            "Gene": gene_name,
            "Module": module,
            "Domain": domain_type,
            "Location": location,
            "AA Sequence": aa_seq,
            "NT Sequence": nt_seq
        })

    return domain_entries

# Step 4: Main script
def main():
    all_data = []
    bgc_ids = get_pks_bgc_ids_from_file("mibig_json_2.0.json")
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False, slow_mo=50)
        page = browser.new_page()
        for bgc_id in bgc_ids:
            all_data.extend(process_bgc(page, bgc_id))
        browser.close()

    if all_data:
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=all_data[0].keys())
            writer.writeheader()
            writer.writerows(all_data)
        print(f"\n✅ Extracted {len(all_data)} ACP/KS entries to {output_file}")
    else:
        print("⚠️ No ACP/KS domains found.")

if __name__ == "__main__":
    main()
