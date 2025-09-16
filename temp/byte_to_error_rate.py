BYTE_TO_ERROR_RATE = []
for i in range(256):
    if i < 33:
        error_rate = 0.0  # Very high error rate for non-printable characters
    else:
        error_rate = 10 ** (-(i - 33) / 10)
    BYTE_TO_ERROR_RATE.append(error_rate)

print(BYTE_TO_ERROR_RATE)
    