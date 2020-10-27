## Hint

### bcftools filter

```
bcftools filter -S . -i 'QUAL>=999' -o filtered.vcf result.vcf
bcftools filter -S . -i 'TYPE=snp' -o filtered.vcf result.vcf
bcftools filter -S . -i 'COUNT(ALT)<=1' -o filtered.vcf result.vcf
bcftools filter -S . -i 'INFO/AD>=10' -o filtered.vcf result.vcf
bcftools filter -S . -i 'QUAL>=999 && TYPE=snp && COUNT(ALT)<=1 && INFO/AD>=10' -o filtered.vcf result.vcf
```

